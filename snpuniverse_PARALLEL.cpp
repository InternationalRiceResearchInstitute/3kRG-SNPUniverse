/*Author: Roven Rommel B. Fuentes
 *TT-Chang Genetic Resources Center, International Rice Research Institute
 *Last Modified: Oct 07,2016 (Added comments)
 *
 *SNP_universe with multithreading by samples and integrated alignment of references

 *COMPILE:
	g++ -o snpuniverse snpuniverse2.cpp -lgzstream -L/home/rfuentes/gzstream -lz -std=c++0x -lpthread
	./snpuniverse -p /shared/data/3kgenomes/newvcf/nipponbare/ -x /data/rfuentes/Nipponbare/3kGenomes_NB -a /home/rfuentes/scaffold_lengths/IRGSP-1.0_genome.out -b /home/rfuentes/feature_ids/nipponbare.dsv -i > /data/rfuentes/Nipponbare/log.txt
 *NOTE: Check NUMTHREADS before running this program.
 *QUAL==30 filter is applied
 *This code assumes that the VCF is emit-all-sites, not needing to enumerate bases in REF (if REF.size()>1 && ALT=.)

 *TO_UPDATE: Migrate all SNPs from INDELuniverse to SNPuniverse
 *         : Current version requires some portion(reduction for multiple REF) of the code to be commented out when running NB. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sys/resource.h>
#include <unordered_map>
#include <vector>
#include <gzstream.h>
#include <time.h>
#include <pthread.h>
#include <sched.h> //for CPU affinity
#include <unistd.h> //for sysconf

using namespace std;

struct chromData{
    int size;
    char* snp;
    int* count;
    char* ref;
};

vector<string> vcf_list;
vector<string> chrom;
unordered_map<string,int> samples;
unordered_map<string,int> chr;
unordered_map<string,int> featureid;
char **refgeno;

struct threadData{
    string dirpath;
    chromData cdata[12730]; //max number of scaffolds
    double* depth;
    double* qs;
    int tid;
};

static const char *options="p:P:x:X:U:u:A:a:B:b:C:c:f:F:mMiI";
static string outfile;
static string path;
static string uniqlist_path;
static string scaffold_size;
static string feature_id;
static string scaffold_algnmnt;
int SIZE=45000000;
int ARRAYSIZE=20000000;
int NUMTHREADS = 25;
int DEPTH = 1;
int COUNTALLSNP=0;
int maxscaf=0;
double offset_snpid=0;
double scafoffset[12730];
bool file_pref=false;
bool indeluni=false;
int last_c,last_p; //last chr+pos to include last indel pos(mono)

void parseArgs(int argc, char**argv){
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
	case 'p':
	case 'P': path = optarg; break; //complete path to the input
	case 'x':
	case 'X': outfile = optarg; break; //output file
        case 'U': 
        case 'u': uniqlist_path = optarg; break; //list of positions from other runs; this parameter is used when upgrading or merging SNP_universes
        case 'A': 
	case 'a': scaffold_size = optarg; break; //list of scaffold sizes 
	case 'B': 
	case 'b': feature_id = optarg; break; //list of scaffold IDs
        case 'C': 
	case 'c': scaffold_algnmnt = optarg; break; //list of scaffold alignment against other reference/s
	case 'f': 
	case 'F': offset_snpid = atof(optarg); break; //sum of lengths of all prev genomes
        case 'm':
        case 'M': file_pref = true; break; //if VCF file names has no prefix
	case 'i':
        case 'I': indeluni = true; break; //get indel universe
	default: break;
        } // switch
    } // while
} // parseArgs

char getcall(char newval, char oldval){
    if(newval=='A'){ //new SNP
	if(oldval=='B') return 'C'; //SNP + Indel
	else if(oldval=='C') return oldval; //retain
	else if(oldval=='X') return oldval; //aligned to previous other reference
	else return 'A'; //SNP
    }else if(newval=='B'){ //newval=='B' ;  new Indel
	if(oldval=='A') return 'C'; //SNP + Indel
	else if(oldval=='C') return oldval; 
	else if(oldval=='X') return oldval; //aligned to previous other reference
	else return 'B';
    }else if(newval=='C'){
	if(oldval=='X') return oldval; //always ingore pos with 'X'
	else return newval;
    }else if(newval=='X'){
	return newval;
    }else{ //0 or initial value; ignore
	return oldval; //retain old call
    }
}

int checkAlt(string ref,string alt,void *thread_data,int pos,int chrpos){
    //return value(max deletion) is needed to reject the GATK emitted positions that intersect with indels
    //this function does not consider the actual genotype in the union
    threadData *t = (threadData*) thread_data;
    int temp=0; 

    if(!alt.compare(".")){ return pos;} //monomorphic
    else if(alt.size()==1 && alt[0]!=ref[0]){ //alt[0]!=ref[0]){ //SNP
	t->cdata[chrpos].snp[pos]=getcall('A',t->cdata[chrpos].snp[pos]);
	t->cdata[chrpos].ref[pos]=ref[0];
	if(!indeluni) t->cdata[chrpos].count[pos]++;
	return pos;
    }

    if(ref.size()>1){ //long ref-> indel, except when ALT='.'
	t->cdata[chrpos].snp[pos]=getcall('B',t->cdata[chrpos].snp[pos]); 
	if(indeluni) t->cdata[chrpos].count[pos]++;
	return (pos+ref.size()-1); //return the longest deletion(always equal to REF.size()-1)
    }else if(alt.size()>1){ //single-base REF
	temp=alt.find_first_of(",",temp+1); //check for multiple ALTs
        if(temp==alt.npos){ //INSERTION
	    t->cdata[chrpos].snp[pos]=getcall('B',t->cdata[chrpos].snp[pos]); 
	    if(indeluni) t->cdata[chrpos].count[pos]++;
	}else{ //Multi-ALT
	    if(alt.size()==3){
		t->cdata[chrpos].snp[pos]=getcall('A',t->cdata[chrpos].snp[pos]); //SNPs
		t->cdata[chrpos].ref[pos]=ref[0];
		if(!indeluni) t->cdata[chrpos].count[pos]++;
	    }else{
		t->cdata[chrpos].snp[pos]=getcall('B',t->cdata[chrpos].snp[pos]); //Insertions
		if(indeluni) t->cdata[chrpos].count[pos]++;
	    }
	}
    }
    return pos; //multiple ALTs
}

int intersecIndel(string ref,string alt,int pos){
    int temp;
    if(!alt.compare(".")){ return pos;} //monomorphic and SNP

    if(ref.size()>1){ //long ref-> indel, except when ALT='.'
	return (pos+ref.size()-1); //return the longest deletion(always equal to REF.size()-1)
    }
    return pos;
}

int locateSNP_1(string filepath,void *thread_data){
    return 0;
}

int locateSNP_2(string filepath,void *thread_data){
    threadData *t = (threadData*) thread_data;
    string linestream,temp,alt,qualstr;
    int idx1=0,idx2=0,pos=0,chrpos=0,ignore=0,prev=0,oldchr=0;
    float qual=0;
    string ref;
    
    igzstream fp(filepath.c_str());
    if (!fp.good()) { //check input file
	printf("ERROR: Failed to open the input file %s", filepath.c_str());
	return 1;
    }
    
    for(int x=0; std::getline(fp,linestream);x++)
    {
      	if(linestream[0]!='#'){
	    idx1 = linestream.find_first_of("\t",idx1+1); //first column
            temp=linestream.substr(0,idx1); //get the chrom number
            idx2 = linestream.find_first_of("\t",++idx1); //second column
            pos = atoi(linestream.substr(idx1,idx2-idx1).c_str()); //get the SNP pos
	   
            if(chr.find(temp)!=chr.end()) chrpos = chr.find(temp)->second; 
	    else { //excluse missing chrom
		//printf("Scaffold not found %s %d.\n",temp.c_str(),pos); 
		//prev=pos; 
		idx1=0; 
		continue;
	    }

	    //Check if genotype = "./." or snp[pos]=='X'
	    if(!strcmp(linestream.substr(linestream.size()-3,3).c_str(),"./.") || t->cdata[chrpos].snp[pos]=='X'){ 
		idx1=0;
		continue;		
	    }
	    
	    idx1 = linestream.find_first_of("\t",idx2+1); //skip id column
    	    idx2 = linestream.find_first_of("\t",idx1+1); //ref
    
	    ref=linestream.substr(idx1+1,idx2-idx1-1); 
	    idx1 = linestream.find_first_of("\t",++idx2);
	    alt=linestream.substr(idx2,idx1-idx2); 

	    qualstr = linestream.substr(idx1+1,linestream.find_first_of("\t",idx1+1)-idx1-1);
	    if(qualstr.compare(".") && qualstr.find_first_not_of("0123456789.") == string::npos){ 
			qual=atof(qualstr.c_str());//qual
            }else{
		idx1=0;
		continue;
	    }	

	    if(qual<30){ //QUAL filtering
		idx1=0;
		continue;
	    }
	
	    if(chrpos!=oldchr){ //needs this to re-initialize ignore and prev for new chr/scaffold
		ignore=0;
		prev=0;
	    } 
            oldchr=chrpos; 

	    //if(chrpos==11 && pos==27530704) printf("%d %s\n",x,linestream.c_str());
	    //get genotype!!!!!!!!!!!!!!!!
            if(pos>ignore || (prev==pos && pos==ignore)){ //ignore positions intersecting with deletion 
		ignore=checkAlt(ref,alt,t,pos,chrpos); 
	    } 
	    prev=pos;
    	}
	idx1=0;
    }
     
    fp.close();
    return 0;
}

void *unionVariant(void *thread_data){
    threadData *t = (threadData*) thread_data;
    int temp,maxfiles,filecount;
    filecount = vcf_list.size();
    string tempfile;
    
    if(filecount%NUMTHREADS){
	if(NUMTHREADS<2) {printf("More processors.\n"); exit(EXIT_FAILURE);}
	maxfiles = filecount/NUMTHREADS;
	temp=t->tid*(maxfiles);
	if(t->tid==(NUMTHREADS-1)) maxfiles+=filecount%NUMTHREADS; //adds the remainder to the last processor
    }else{
	maxfiles = filecount/NUMTHREADS;
 	temp=t->tid*(maxfiles);
    }
    printf("thread: %d\tmaxfile:%d\n",t->tid,maxfiles);
    for(int i=0;i<maxfiles;i++){
	//tempfile = "gzip -t " + vcf_list[i+temp] + " 2>&1 >> log.txt";
	//system(tempfile.c_str());
	locateSNP_2(vcf_list[i+temp],thread_data);
	printf("%d %s\n",t->tid,vcf_list[i+temp].c_str());
	fflush(stdout);
    }
			
}

void readFolder(string dirpath){
    DIR *dp;
    string temp,exten,filename,curpath;
    struct dirent *ep;
    int dotpos;
    
    dp = opendir(dirpath.c_str());
    if(dirpath[dirpath.size()-1]!='/') curpath = dirpath + '/';
    else curpath = dirpath;
    struct stat filestat;
    if(dp!=NULL){
	while(ep = readdir(dp)){
	    if(!strcmp(ep->d_name,".") || !strcmp(ep->d_name,"..")) continue;
	    temp = curpath + ep->d_name;
            stat(temp.c_str(), &filestat);
	    filename=ep->d_name;
	    if(S_ISDIR(filestat.st_mode)){ //recursively read directories
                dirpath=temp;
		readFolder(dirpath);
                dirpath=curpath;
	    }else if(S_ISREG(filestat.st_mode) || S_ISLNK(filestat.st_mode)){
		//check if vcf
		dotpos = filename.size()-7; 
		if(dotpos>0){
		    exten=filename.substr(dotpos);
		    if(!exten.compare(".vcf.gz")){
			vcf_list.push_back(temp);
			printf("%s\n",temp.c_str());
	            }
		}

	    }
	    temp.clear();
	    filename.clear();
	}
	closedir (dp);
    }else{
	printf("ERROR: Can't find the directory.");
    }
 
}

void *multiprint_2(void *thread_data){
    threadData *t = (threadData*) thread_data;
    string linestream,refstr,alt,temp1,formatfield,formatval,temp2,samname,geno,tempLine; 
    char ref,*token=NULL,tok_ar[40],buffer[15];
    int idx1=0,idx2=0,idx3=0,idx4=0,pos=0,chrpos=0,DP=0,filecount=samples.size();
    int AD[128],GID=0,prev=0,commaidx=0;
    char alleles[4];
    float qual;
    double pos_id=0;//initialize to size of previous reference/s

    int first=0,last=0,maxfiles;
    if(filecount%NUMTHREADS){
	if(NUMTHREADS<2) {printf("More processors.\n"); exit(EXIT_FAILURE);}
	maxfiles = filecount/NUMTHREADS;
	first=t->tid*(maxfiles);
	if(t->tid==(NUMTHREADS-1)) maxfiles+=filecount%NUMTHREADS;
    }else{
	maxfiles = filecount/NUMTHREADS;
 	first=t->tid*(maxfiles);
    }

    last = first+maxfiles-1; //index of the last sample in the subset
    //printf("thread:%d %d %d\n",t->tid,first, last);
    FILE *output1;
    if(indeluni) temp1=outfile+"_INDEL_GT_"+ to_string(static_cast<long long>(t->tid+1)) +".txt"; 
    else temp1=outfile+"_SNP_GT_"+ to_string(static_cast<long long>(t->tid+1)) +".txt"; 
    output1=fopen(temp1.c_str(),"w"); 
    
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(t->tid+1,&cpuset);

    sched_setaffinity(0,sizeof(cpuset),&cpuset);

    //PRINT Sample Names with GID/BoxCode
    for(int i=0;i<filecount;i++){ 
	idx1=0;
	if(file_pref){
	    idx1=vcf_list[i].find_last_of("/")+1;
	    samname = vcf_list[i].substr(idx1,vcf_list[i].size()-idx1-7);
	}else{
	    idx1=vcf_list[i].find_first_of("_",idx1); //ignore 1st "_"
	    idx1=vcf_list[i].find_first_of("_",idx1+1)+1; //ignore 2nd "_"
	    samname = vcf_list[i].substr(idx1,vcf_list[i].size()-idx1-11);
	}
        
        //sample list
	if(samples.find(samname)!=samples.end()) 
	    fprintf(output1,">%d %s\n",samples.find(samname)->second,samname.c_str()); 
    }
    if(indeluni) fprintf(output1,"SampleID\tINDEL_ID\tPos\tRef\tAlt\tQS\tGenotype\n");
    else fprintf(output1,"SampleID\tSNP_ID\tChr\tPos\tA\tC\tG\tT\tQS\n");//header
	
    
    int snpcount=0,all1=0,all2=0, ignore=0,oldchr=first;
    bool aVariant = false;
    
    int puncidx1,puncidx2;
    for(int i=first;i<=last;i++){
        igzstream fp(vcf_list[i].c_str());
    	if (!fp.good()) { //check input file
	    printf("ERROR: Failed to open the input file %s", vcf_list[i].c_str());
	    exit(EXIT_FAILURE);
    	}
	
        
  	//get sample name
	idx1=0;
	if(file_pref){
	    idx1=vcf_list[i].find_last_of("/")+1;
	    samname = vcf_list[i].substr(idx1,vcf_list[i].size()-idx1-7);
	}else{
	    idx1=vcf_list[i].find_first_of("_",idx1); //ignore 1st "_"
	    idx1=vcf_list[i].find_first_of("_",idx1+1)+1; //ignore 2nd "_"
	    samname = vcf_list[i].substr(idx1,vcf_list[i].size()-idx1-11);
	}
        
        if(samples.find(samname)!=samples.end()) GID = samples.find(samname)->second;
	else{ printf("The sample has no GID:%s\n",samname.c_str()); continue;};
      
        printf("Printing:\t%s\n",samname.c_str());
	fflush(stdout);
        
        idx1=0;
        ignore=0;
        prev=0;
	pos=0;

    	for(int x=0;getline(fp,linestream);x++){
	    if(linestream[0]!='#'){ 
		idx1=DP=0;
		prev=pos;
                idx1 = linestream.find_first_of("\t",idx1+1); //first column
             	temp1 = linestream.substr(0,idx1); //get the chrom number
            	idx2 = linestream.find_first_of("\t",++idx1); //second column
            	pos = atoi(linestream.substr(idx1,idx2-idx1).c_str()); //get the SNP pos
                if( chr.find(temp1)!=chr.end()) chrpos = chr.find(temp1)->second;
		else{ //excluse missing chrom
		    continue;
		}

   	        if(chrpos!=oldchr){ //needs this to re-initialize ignore and prev for new chr/scaffold
		    ignore=0;
		    prev=0;
		    pos=0;
		} 
                oldchr=chrpos; 

		//Check if genotype = "./." or snp[pos]=='X'
		if(!strcmp(linestream.substr(linestream.size()-3,3).c_str(),"./.") || t->cdata[chrpos].snp[pos]=='X') continue;		
		
                idx1 = linestream.find_first_of("\t",++idx2)+1; //skip ID
            	idx2 = linestream.find_first_of("\t",idx1); //ref
		//get REF
		refstr=linestream.substr(idx1,idx2-idx1);
		//get ALT
		idx1=linestream.find_first_of("\t",++idx2);
	 	alt=linestream.substr(idx2,idx1-idx2); //alt

                //get QUAL
		idx2=linestream.find_first_of("\t",++idx1);
		temp1 = linestream.substr(idx1,idx2-idx1);
		if(idx2!=idx1 && temp1.find_first_not_of("0123456789.") == string::npos){ 
		    qual=atof(temp1.c_str());//qual
		}else{
		    qual=0;
		    //printf("WARNING: Invalid QUAL for \"%s\" at position: %.0f:\"%s\"\n",samname.c_str(),scafoffset[chrpos]+pos,temp1.c_str()); 
		    continue;
		}

            	if(t->cdata[chrpos].snp[pos]=='A' && !indeluni){ //COUNTALLSNP++;
		    //QUAL filtering: remove low QUAL pos
                    if(qual<30) continue;

                    //CHECK for intersecting INDEL
                    if(pos>ignore || (prev==pos && pos==ignore)){
                    	ignore=intersecIndel(refstr,alt,pos);
                    	//ignore positions intersecting with deletion
                    	if(ignore>pos) continue; //skip next pos until last base of the deletion
                    }else{ //ignore positions intersecting with deletion
                    	continue;
                    }

		    //get REF
		    if(refstr.size()==1){ 
		        ref=refstr[0];
		        alleles[0]=refstr[0];
		    }else{
			if(alt[0]=='.'){ //REF.size()>1 and ALT='.'
			    ref=refstr[0];
			}else{  
			   printf("ERROR: Indel is mixed. %d %s\n",GID,linestream.c_str());
			   aVariant=false;
			   continue;
			}
		    }

   
		    //get ALT
		    alleles[1]=alt[0];
		    if(alt.size()>1) alleles[2]=alt[2]; 

		    idx1=idx2;
	    	    for(int y=0;y<2;y++)idx1 = linestream.find_first_of("\t",idx1+1); //3columns
		    idx2 = linestream.find_first_of("\t",++idx1);
		    //format field IDs
    	    	    formatfield=linestream.substr(idx1,idx2-idx1); 
		    //format value 
		    formatval=linestream.substr(idx2+1); 

		    fprintf(output1,"%d\t%.0f\t%s\t%d\t",GID,scafoffset[chrpos]+pos,chrom[chrpos].c_str(),pos-1); //PRINT variety and scaffold/chrom ID
 		    //print allele depth
		    idx1=idx3=0;
	            while(idx1!=formatfield.npos){   
			idx2 = formatfield.find_first_of(":\0",++idx1); 
			temp2=formatfield.substr(idx1,idx2-idx1); //get field ID
			idx4 = formatval.find_first_of(":\0",++idx3); 
			if(idx1==1){
				//get genotype
				all1=formatval[0]-'0'; 
				all2=formatval[2]-'0';
			}else if(!temp2.compare("AD")){
                                AD['A']=AD['T']=AD['C']=AD['G']=0;
				string ADval = formatval.substr(idx3,idx4-idx3+1);
				puncidx1 = ADval.find_first_of(",:");
				AD[ref]=atoi(ADval.substr(0,puncidx1).c_str());
				if(puncidx1!=ADval.size()-1){	
				    puncidx2 = ADval.find_first_of(",:",++puncidx1);
				    AD[alt[0]]=atoi(ADval.substr(puncidx1,puncidx2-puncidx1).c_str());
				    if(puncidx2!=ADval.size()-1){
					puncidx1 = ADval.find_first_of(",:",++puncidx2);
				        AD[alt[2]]=atoi(ADval.substr(puncidx2,puncidx1-puncidx2).c_str());
				    }	
				}else{
				    printf("ERROR: Invalid AD/DP %d %s.\n",GID,linestream.c_str());
				}
				aVariant = true;
				fprintf(output1,"%d\t%d\t%d\t%d\t%.2f\t%c\t%c\n",AD['A'],AD['C'],AD['G'],AD['T'],qual,alleles[all1],alleles[all2]);
						       
			}else if(!temp2.compare("DP")){
			    temp1 = formatval.substr(idx3,idx4-idx3);
			    if(idx4!=idx3 && temp1.find_first_not_of("0123456789") == string::npos) 
				DP=atoi(temp1.c_str());
			    else{
			        DP=0;
				printf("WARNING: Invalid depth(DP) for \"%s\" at position: %.0f:\"%s\"\n",samname.c_str(),scafoffset[chrpos]+pos,temp1.c_str()); 
			    }
			}
			idx1=idx2;
			idx3=idx4;
            	    }
		    if(!aVariant){ //Not a SNP but the position is included in the union
			AD['A']=AD['T']=AD['C']=AD['G']=0;
			AD[ref]=DP; 
			//fprintf(output1,"%d\t%d\t%d\t%d\t%.2f\t%c\n",AD['A'],AD['C'],AD['G'],AD['T'],qual,alleles[all1]);
			fprintf(output1,"%d\t%d\t%d\t%d\t%.2f\t%c\t%c\n",AD['A'],AD['C'],AD['G'],AD['T'],qual,alleles[all1],alleles[all2]);
		    } 
		    if(DP<ARRAYSIZE){
		    	if(DP>0) t->depth[(int)DP]++;
                    	else t->depth[0]++;
		    }else printf("WARNING:Bin for DP is overflowing in pos %d!\n",pos);

		    if(qual<ARRAYSIZE){
		    	if(qual>0) t->qs[(int)qual]++;
                    	else t->qs[0]++;
		    }else printf("WARNING:Bin for QUAL is overflowing in pos %d!\n",pos);
	       	}else if(indeluni && (t->cdata[chrpos].snp[pos]=='B' || t->cdata[chrpos].snp[pos]=='C')){
		    if(!tempLine.empty() && prev!=pos){
			fprintf(output1,"%s\n",tempLine.c_str());
		    }
		    tempLine.clear();

                    //QUAL filtering: remove low QUAL pos
                    if(qual<30) continue;

                    //CHECK for intersecting INDEL
                    if(pos>ignore || (prev==pos && pos==ignore)){
                        ignore=intersecIndel(refstr,alt,pos); //ignore pos after the DEL anchor 
                    }else{ //ignore positions intersecting with deletion
                        continue;
                    }

               	    idx1=idx2;
	    	    for(int y=0;y<2;y++)idx1 = linestream.find_first_of("\t",idx1+1); //3columns
		    idx2 = linestream.find_first_of("\t",++idx1);
		    //format field IDs
    	    	    formatfield=linestream.substr(idx1,idx2-idx1); 
		    //format value 
		    formatval=linestream.substr(idx2+1);

		    geno = formatval.substr(0,formatval.find_first_of(":\n",0)); //genotype

		    if(alt[0]=='.'){  //This is necessary to skip anchor pos and only print the variant
			sprintf(buffer,"%.2f",qual);
			if(refstr.size()>1) refstr[1]='\0'; //concatenated monomorphic genotypes
			tempLine = to_string(static_cast<long long>(GID)) +"\t"+ to_string(static_cast<long long>(scafoffset[chrpos]+pos)) +"\t"+ chrom[chrpos].c_str();
			tempLine += "\t"+ to_string(static_cast<long long>(pos-1)) +"\t"+ refstr.c_str() +"\t"+ alt.c_str() +"\t"+ buffer +"\t"+ geno.c_str();
		    }else{
			fprintf(output1,"%d\t%.0f\t%s\t%d\t%s\t%s\t%.2f\t%s\n",GID,scafoffset[chrpos]+pos,chrom[chrpos].c_str(),pos-1,refstr.c_str(),alt.c_str(),qual,geno.c_str());
		    }
		    
		    if(qual<ARRAYSIZE){
		    	if(qual>0) t->qs[(int)qual]++;
                    	else t->qs[0]++;
		    }else printf("WARNING:Bin for QUAL is overflowing in pos %d!\n",pos);
		
	        }

		if( last_p==pos-1 && last_c==chrpos && !tempLine.empty()){
	 	    fprintf(output1,"%s\n",tempLine.c_str());	
		    tempLine.clear();
		}
		
		aVariant=false;
	    }
    	} //end of loop
        //printf("ALL SNP positions in %s: %d\n",samname.c_str(),COUNTALLSNP);
        //COUNTALLSNP=0;
	fp.close();
	 printf("Done printing:\t%s\n",samname.c_str());
    }   
    fclose(output1);
}

int main(int argc, char **argv)
{
    const rlim_t STACK_SIZE = 1000*1024*1024*4lu; 
    struct rlimit rl;
    rl.rlim_cur = STACK_SIZE;
    int ret = setrlimit(RLIMIT_STACK,&rl);
    double *tally1;
    time_t start, end;
    pthread_t  threads[NUMTHREADS];
    unordered_map<string,int> boxcode_map,iris_map; 

    time(&start);
    parseArgs(argc,argv);
    if (path.empty() || outfile.empty() || scaffold_size.empty()) {
        cerr << "Usage:\n" << *argv
                  << " -p path to directory\n"
		  << " -i chromosome index\n"
		  << " -x output file\n"
		  << " -c compressed VCFs (vcf.gz)\n"
		  << " -u uniq list for universe merging\n"
	    	  << " -a scaffold names and sizes\n"
		  << " -b scaffold IDs\n"
		  << " -c scaffold alignment against other refs\n"
                  << endl;
        return 1;
    }
    
    //load list of GIDs
    int pos1,pos2,pos3,gid;
    int idx1=0,idx2=0,cntr=0;
    threadData thread_data[NUMTHREADS]; 
    int *sizes = (int*)calloc(12730,sizeof(int)); 
    string samname,linestream,gidstr,errpath;
    
    if(path[path.size()-1]=='/') path = path.substr(0,path.size()-1); //remove the last '/' from the path
    ifstream list("/home/rfuentes//SNPuniverse/3k_ID_BOXCODE_IRIS.txt");
    if(!list.is_open()){
	printf("ERROR: Cannot open the GID list.");
	return 1;
    }
    
    for(int x=0;getline(list,linestream)!=NULL;x++){
	gid=0;
        pos1 = linestream.find_first_of(",");
        pos2 = linestream.find_first_of(",",pos1+1);
	gidstr = linestream.substr(pos1+1,pos2-pos1-1);
	if(gidstr.find_first_not_of("0123456789.") == string::npos) gid = atoi(gidstr.c_str());
	boxcode_map.insert(pair<string,int>(linestream.substr(0,pos1),gid));
	//printf("%s ",linestream.substr(0,pos1).c_str());
	pos1 = linestream.find_first_of(",",pos2+1); //save IRIS_name
	iris_map.insert(pair<string,int>(linestream.substr(pos2+1,pos1-pos2-1),gid));
	//printf("%d %s\n",gid,linestream.substr(pos2+1,pos1-pos2-1).c_str());
    }
        
    ifstream scafsize(scaffold_size.c_str());
    if (!scafsize.is_open()) { //check input file
	printf("ERROR: Can't open list of scaffold size.");
	exit(EXIT_FAILURE);
    }
   
       
    /*ifstream alignmnt(scaffold_algnmnt.c_str());
    if (!alignmnt.is_open() && !scaffold_algnmnt.empty()) { //check input file
	printf("ERROR: Can't open alignment file.");
	exit(EXIT_FAILURE);
    }*/

    ifstream feat_id(feature_id.c_str());
    if (!feat_id.is_open()) { //check input file
	printf("ERROR: Can't open feature ID list.");
	exit(EXIT_FAILURE);
    }
    
    //read chrom names and sizes
    int genomesize=0;
    scafoffset[0]=0;
    for(cntr=0;getline(scafsize,linestream);cntr++){
	idx1 = linestream.find_first_of("\t");
	chr[linestream.substr(0,idx1)] = cntr; //read all chrom/scaffold
	chrom.push_back(linestream.substr(0,idx1)); 
	sizes[cntr] = atoi(linestream.substr(idx1).c_str());
	genomesize += sizes[cntr];
	if(cntr>0) scafoffset[cntr]=scafoffset[cntr-1]+sizes[cntr-1];
	//printf("%s %d\n",linestream.substr(0,idx1).c_str(),sizes[cntr]);	
    }
    maxscaf = cntr; 
    for(int i=0;i<maxscaf;i++){
	scafoffset[i] += offset_snpid; //adding size of prev genomes
    }
    printf("TOTAL genome size: %d.\n",genomesize);

    //read feature_ids
    getline(feat_id,linestream); //remove header
    for(int i=0;getline(feat_id,linestream);i++){
	featureid[linestream.substr(linestream.find_last_of("\t")+1)]=atoi(linestream.substr(0,linestream.find_first_of("\t")).c_str());
	//printf("%s %d\n",linestream.substr(linestream.find_last_of("\t")+1).c_str(),featureid[linestream.substr(linestream.find_last_of("\t")+1)]);
    }
    
    string dirpath=path;
    readFolder(dirpath);
    
    tally1=(double*)calloc(vcf_list.size()+1,sizeof(double));
     
    /*THREADING of the comparison*/
    string tempname;
    int first,last,scafidx;	
    for(int i=0;i<NUMTHREADS; i++){
        thread_data[i].dirpath=path;
        for(int j=0;j<maxscaf;j++){
	    thread_data[i].cdata[j].size = sizes[j];
	    thread_data[i].cdata[j].snp = (char*)calloc(sizes[j],sizeof(char));
	    if(thread_data[i].cdata[j].snp==NULL){ printf("Not enough space for snp array."); return 1;}
	    thread_data[i].cdata[j].count = (int*)calloc(sizes[j],sizeof(int));
	    if(thread_data[i].cdata[j].count==NULL){ printf("Not enough space for count array."); return 1;}
	    thread_data[i].cdata[j].ref = (char*)calloc(sizes[j],sizeof(char));
	    if(thread_data[i].cdata[j].ref==NULL){ printf("Not enough space for ref array."); return 1;}
	}

        //read alignment and flag positions that intersect with other reference/s
	/*if(!scaffold_algnmnt.empty()){
	    for(int j=0;getline(alignmnt,linestream);j++){
	    	idx1=linestream.find_first_of("\t")+1;
	    	idx2=linestream.find_first_of("\t",idx1);
	    	tempname=linestream.substr(idx1,idx2-idx1);
	    	idx1=linestream.find_first_of("\t",++idx2);
	    	first =atoi(linestream.substr(idx2,idx1-idx2).c_str());
	    	last =atoi(linestream.substr(++idx1).c_str());
	    	scafidx = chr.find(tempname)->second; //find idx for scaffold array
	    	//if(j==0)printf("%d %d\n",first,last);
	    	for(int k=first;k<=last;k++){
		    thread_data[i].cdata[scafidx].snp[k] = 'X'; 
	            //printf("%c\n",thread_data[i].cdata[scafidx].snp[k]);
	     	}
	    	//printf("%s %d %d\n",tempname.c_str(),first,last);
	    }
	    alignmnt.clear();
	    alignmnt.seekg(0,ios::beg);
	}*/
    	
	thread_data[i].depth = (double*)calloc(ARRAYSIZE,sizeof(double));
        thread_data[i].qs = (double*)calloc(ARRAYSIZE,sizeof(double));
	
	if(thread_data[i].depth==NULL){ printf("Not enough space for depth array."); return 1;}
        if(thread_data[i].qs==NULL){ printf("Not enough space for qs array."); return 1;}
    	thread_data[i].tid = i;
      
	pthread_create(&threads[i],NULL,unionVariant, (void*) &thread_data[i]); 
    }
   
    for(int i=0;i<NUMTHREADS;i++){
	pthread_join(threads[i],NULL);
    }
    
    //ADDING VCFs from other directory to the set
    /*for(int i=0;i<NUMTHREADS;i++){ 
	thread_data[i].dirpath="/storage2/halcon3k/bam_vcf";
	pthread_create(&threads[i],NULL,readFolder, (void*) &thread_data[i]); 
    }
  
    for(int i=0;i<NUMTHREADS;i++){
	pthread_join(threads[i],NULL);
    }*/

    /*if(!uniqlist_path.empty()){ 
        //MERGING of UNIVERSES
    	for(int i=0;i<NUMTHREADS; i++){
            pthread_create(&threads[i],NULL,readUniqList, (void*) &thread_data[i]);
   	}

    	for(int i=0;i<NUMTHREADS;i++){
            pthread_join(threads[i],NULL);
    	}
    }*/   
    
    //union of all thread_data
    for(int i=1;i<NUMTHREADS;i++){
	for(int j=0;j<maxscaf;j++){
	    for(int k=0;k<sizes[j];k++){ 
	 	//save union all scaffold union to thread_data[0]
	    	thread_data[0].cdata[j].snp[k]=getcall(thread_data[i].cdata[j].snp[k],thread_data[0].cdata[j].snp[k]);
		thread_data[0].cdata[j].count[k]+=thread_data[i].cdata[j].count[k];
		if(thread_data[i].cdata[j].ref[k]!=0){
		    thread_data[0].cdata[j].ref[k] = thread_data[i].cdata[j].ref[k];
		}
	    }
	}
    }
    //update thread_data
    for(int i=1;i<NUMTHREADS;i++){
	if(i==NUMTHREADS-1){ // get location of last variant pos
	    for(int j=0;j<maxscaf;j++){
	    	for(int k=0;k<sizes[j];k++){ 
	 	    //save union all scaffold union to thread_data[0]
	    	    thread_data[i].cdata[j].snp[k]=thread_data[0].cdata[j].snp[k];
		    if(thread_data[i].cdata[j].snp[k]=='B' || thread_data[i].cdata[j].snp[k]=='C'){
		        last_c=j;
			last_p=k;
		    }
	        }
	    }
	}else{
	    for(int j=0;j<maxscaf;j++){
	    	for(int k=0;k<sizes[j];k++){ 
	 	    //save union all scaffold union to thread_data[0]
	    	    thread_data[i].cdata[j].snp[k]=thread_data[0].cdata[j].snp[k];
	        }
	    }
	}
    }
    
    int setsize=vcf_list.size();
    for(int i=0;i<setsize;i++){ 
	//GET IDs for each VCF/sample
	idx1=0; 
	if(file_pref){
	    idx1=vcf_list[i].find_last_of("/")+1;
	    samname = vcf_list[i].substr(idx1,vcf_list[i].size()-idx1-7);  
	}else{
  	    idx1=vcf_list[i].find_first_of("_",idx1); //ignore 1st "_"
	    idx1=vcf_list[i].find_first_of("_",idx1+1)+1; //ignore 2nd "_"
	    samname = vcf_list[i].substr(idx1,vcf_list[i].size()-idx1-11);  
	}
        
        if(samname.find("IRIS")==string::npos){ //no IRIS name
	   if(boxcode_map.find(samname)!=boxcode_map.end())
           	samples.insert(pair<string,int>(samname,boxcode_map.find(samname)->second));
	}else{ //with IRIS name
	   if(iris_map.find(samname)!=iris_map.end())
	   	samples.insert(pair<string,int>(samname,iris_map.find(samname)->second));
	}
    }
 
    time(&end);
    boxcode_map.clear();
    iris_map.clear();
    printf("Union Time: %.f sec\n",difftime(end,start));
    fflush(stdout); 
    
    //PRINTING universe
    for(int i=0;i<NUMTHREADS; i++){
	pthread_create(&threads[i],NULL,multiprint_2, (void*) &thread_data[i]);
    }

    for(int i=0;i<NUMTHREADS;i++){
	pthread_join(threads[i],NULL);
    } //END Printing of universe
    
    //PRINTING Statistics
    FILE *output1,*output2,*output3,*output4,*output5;
    string TYPE;
    if(indeluni) TYPE="INDEL";
    else TYPE="SNP";
    string temp=outfile+"_"+TYPE+"_POS_vs_VARcount.txt"; 
    output1=fopen(temp.c_str(),"w");  
    temp=outfile+"_"+TYPE+"_POS_vs_DEPTH.txt"; 
    output2=fopen(temp.c_str(),"w"); 
    temp=outfile+"_"+TYPE+"_POS_vs_QS.txt"; 
    output3=fopen(temp.c_str(),"w"); 
    temp=outfile+"_"+TYPE+"_SUPPORT.txt";
    output4=fopen(temp.c_str(),"w");
    temp=outfile+"_"+TYPE+"_POS.txt";
    output5=fopen(temp.c_str(),"w");

    //union of all stats thread_data
    for(int i=1;i<NUMTHREADS;i++){
	for(int j=0;j<ARRAYSIZE;j++){
	    thread_data[0].depth[j]+=thread_data[i].depth[j];
	    thread_data[0].qs[j]+=thread_data[i].qs[j];
	}
    }
    
    /*for(int i=0;i<NUMTHREADS;i++){
	for(int j=0;j<maxscaf;j++){
	    for(int k=0;k<sizes[j];k++){ 
		if(thread_data[i].cdata[j].count[k]) fprintf(output4,"%s\t%d\t%d\n",thread_data[0].chrom[j].c_str(),k,thread_data[i].cdata[j].count[k]);
		tally1[thread_data[i].cdata[j].count[k]]++;
	    	//thread_data[0].cdata[j].count[k]+=thread_data[i].cdata[j].count[k];
	    }
	}
    }*/

    fprintf(output1,"#ofVar\t#ofSNPpos\n");
    fprintf(output2,"Depth\t#ofCalls\n");
    fprintf(output3,"QS\t#ofCalls\n");
    fprintf(output4,"ChrID\tPos\tSupport\n");
    if(indeluni) fprintf(output5,"INDEL_ID\tChr_ID\tChr\t(Pos-1)\n");  
    else fprintf(output5,"SNP_ID\tChr_ID\tChr\t(Pos-1)\tRef\n");  
    int variantcount=0;
    
    if(indeluni){
	for(int j=0;j<maxscaf;j++){
	    for(int k=0;k<sizes[j];k++){ 
	    	if(thread_data[0].cdata[j].snp[k]=='B' || thread_data[0].cdata[j].snp[k]=='C'){
		    variantcount++;
 		    fprintf(output4,"%s\t%d\t%d\n",chrom[j].c_str(),k,thread_data[0].cdata[j].count[k]);
		    fprintf(output5,"%.0f\t%d\t%s\t%d\n",scafoffset[j]+k,featureid[chrom[j]],chrom[j].c_str(),k-1); //pos-1
		    tally1[thread_data[0].cdata[j].count[k]]++; //distribution of supports
	        }
    	    }
            printf("SNP count:\t%s\t%d SNPs.\n",chrom[j].c_str(),variantcount);
	    variantcount=0;
    	}
    }else{
        for(int j=0;j<maxscaf;j++){
	    for(int k=0;k<sizes[j];k++){ 
	    	if(thread_data[0].cdata[j].snp[k]=='A' || thread_data[0].cdata[j].snp[k]=='C'){ //also add SNPs in anchors/ 'C'
		    variantcount++;
 		    fprintf(output4,"%s\t%d\t%d\n",chrom[j].c_str(),k,thread_data[0].cdata[j].count[k]);
		    fprintf(output5,"%.0f\t%d\t%s\t%d\t%c\n",scafoffset[j]+k,featureid[chrom[j]],chrom[j].c_str(),k-1,thread_data[0].cdata[j].ref[k]); //pos-1
		    tally1[thread_data[0].cdata[j].count[k]]++; //distribution of supports
	        }
    	    }
            printf("SNP count:\t%s\t%d SNPs.\n",chrom[j].c_str(),variantcount);
	    variantcount=0;
    	}
    }
    
    for(int i=1;i<=vcf_list.size();i++) 
	if(tally1[i]) fprintf(output1,"%d\t%.0f\n",i,tally1[i]);

    for(int j=1;j<ARRAYSIZE;j++){
	 if(thread_data[0].depth[j]) fprintf(output2,"%d\t%.0f\n",j,thread_data[0].depth[j]);
         if(thread_data[0].qs[j]) fprintf(output3,"%d\t%.0f\n",j,thread_data[0].qs[j]); 
    }

    time(&end);
    printf("Running Time: %.f sec\n",difftime(end,start));
    fflush(stdout);
    //END PRINTING of Stat
   
    
    //DEALLOCATION
    free(tally1);
    chr.clear();
    samples.clear();
    vcf_list.clear();
    fclose(output1);
    fclose(output2);
    fclose(output3);
    fclose(output4);    
    fclose(output5);    
    /*for(int i=0;i<NUMTHREADS;i++){
	for(int j=0;j<maxscaf;j++){
	    free(thread_data[i].cdata[j].snp);
	    free(thread_data[i].cdata[j].count);
	}
	free(thread_data[i].depth);
	free(thread_data[i].qs);
    }*/
    return 0;
}