#pragma once
#include <msclr/marshal_cppstd.h>
using namespace msclr::interop;
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<iostream>
#include<windows.h>
#include<MMSystem.h>
//#include"winmm.lib"
using namespace std;
std::wstring s2ws(const std::string& s)
{
    int len;
    int slength = (int)s.length() + 1;
    len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0); 
    wchar_t* buf = new wchar_t[len];
    MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
    std::wstring r(buf);
    delete[] buf;
    return r;
}
#define pi 3.142857


#define N 5
#define T 85
#define M 32
long double A[N+1][N+1];
long double B[N+1][M+1];
long double Pi[N+1];
long double AVG_A[10][N+1][N+1]={0};
long double AVG_B[10][N+1][M+1]={0};
long double AVG_Pi[10][N+1]={0};
int O[200];

long double Alpha[200][N+1];
long double Beta[200][N+1];
long double Delta[200][N+1];
long double Gamma[200][N+1];
int shi[200][N+1];
long double Eta[200][N+1][N+1];
long double Universe[200*200][13];
int global_i=1;
int cluster[33][200000];
long double codebook[33][13]={0};
int indexes[33]={0};
//long double P_Star;
char Names[100][100];
char ContactNo[100][100];
char deleted[100][100];
FILE *globalFile=fopen("global.txt","w");
int Gflag=1;
int index;
void initial(){
	FILE *fp,*fp1,*fp2,*fp3;
	fp=fopen("Names.txt","r");
	fp1=fopen("ContactNo.txt","r");
	int count=0;
	fp3=fopen("codebook.txt","r");
	for(int i=1;i<=32;i++){
		for(int j=1;j<=12;j++){
			fscanf(fp3,"%Lf",&codebook[i][j]);
		}
	}
	fclose(fp3);
	char temp[100];
	char temp1[100];
	while(fgets(temp,80,fp)!=NULL && fgets(temp1,80,fp1)!=NULL){
		strcpy(Names[count],temp);
		strcpy(ContactNo[count++],temp1);
	}
	index=count;
	for(int i=0;i<index;i++){
		int l1=strlen(Names[i]);
		Names[i][l1-1]='\0';
		l1=strlen(ContactNo[i]);
		ContactNo[i][l1-1]='\0';
	}
	fclose(fp);
	fclose(fp1);
}
//Centroid Calculation

void CentroidCalc(int size){

	
	//Calculating the Centroid and Updation of Code Books
	for(int i=1;i<=size;i++){
		for(int j=1;j<=12;j++){
			codebook[i][j]=0;
		}
	}

	//Centroid is Calculated and stored in the codebook
	for(int i=1;i<=size;i++){
		for(int j=0;j<indexes[i];j++){
			for(int k=1;k<=12;k++){
				codebook[i][k]=codebook[i][k]+Universe[cluster[i][j]][k];
			}
		}
		for(int k=1;k<=12;k++){
			if(indexes[i]!=0){
				codebook[i][k]/=indexes[i];
			}
		}
	}
	
}


//K-Means Function
void kmeans(int size){
	
	while(size<=32){
		long double prev_distortion=INT_MAX;
		
		//Splitting of Code Books
		
		if(size>1){
			long double tempC[33][13];	
			for(int i=1;i<=size/2;i++){
				for(int j=1;j<=12;j++){
					tempC[i][j]=codebook[i][j]*1.03;
				}
			}
			for(int i=size/2+1;i<=size;i++){
				for(int j=1;j<=12;j++){
					tempC[i][j]=codebook[i-(size/2)][j]*0.97;
				}
			}

			for(int i=1;i<=size;i++){
				for(int j=1;j<=12;j++){
						codebook[i][j]=tempC[i][j];
				}
			}
		}
		while(1){
			long double overall_distortion=0;
			long double tokhurasWeights[12]={1.0,3.0,7.0,13.0,19.0,22.0,25.0,33.0,42.0,50.0,56.0,61.0};
			int temp=-1;
			//Making Cluster Sizes to Zero
			for(int i=1;i<=size;i++){
				indexes[i]=0;
			}

			//Dividing the points in the universe into 8 clusters based on distortion
			for(int i=1;i<=global_i;i++){
				long double min=INT_MAX;
				long double sum=0;
				for(int j=1;j<=size;j++){
					sum=0;
					for(int k=1;k<=12;k++){
						long double diff=Universe[i][k]-codebook[j][k];
						sum+=tokhurasWeights[k-1]*diff*diff;
					}
				
					if(min>sum){
						min=sum;
						temp=j;
					}
				}
				overall_distortion+=min;
				cluster[temp][indexes[temp]++]=i;
				
			}
		
			//Calculating Average Distortion
			long double avg_distortion=overall_distortion/global_i;

			//Centroid Calculation
			CentroidCalc(size);

			//Termination of K-Means Algo 
			if(prev_distortion > avg_distortion){
				if(prev_distortion-avg_distortion<=0.0001){
					break;
				}
			}
			else if(avg_distortion -prev_distortion<=0.0001){
				break;
			}
			prev_distortion=avg_distortion;
			
		}
		size*=2;
	}

}


//LBG Function
void LBG(){
	//Initially placing all the points in Universe into the same cluster
	indexes[1]=global_i;
	for(int i=1;i<=global_i;i++){
		cluster[1][i]=i;
	}
	//Making the initial Code Book as the Centroid of the Universe
	CentroidCalc(1);
	//Calling K-means
	kmeans(1);
}









//SteadyFrames will be stored in a file
void CalculateSteadyFrames(char *str){
	
	FILE *fp,*fp1;
	//STEP 1:FINDING DC-SHIFT
	
	
	//Opening a file for calculating the dc shift
	fp=fopen("DCshift.txt","r");
	//Initially first 5 lines are some of the configuration details so we have to skip those by retreving it into this string
	char InitialStrings[80];
	//This is to know or to keep track of number of samples(or)lines we read
	long int cnt=0;
	//ambient noise
	long double Ambient_noise=0;
	//For skipping the first 5 lines
	while(fgets(InitialStrings,80,fp)!=NULL && cnt<5){
		cnt++;
	}
	cnt=0;
	//to store the sample values temporarily in this variable
	int temp=0;
	//after that cool edit will add some headers at the beginning. so we skip first 300 samples
	while(fscanf(fp,"%d",&temp)!=EOF && cnt<300){
		cnt++;
	}
	//now we calculate the DC-shift by taking average of all the samples in this file
	long double DC_Shift=0;
	cnt=0;
	while(fscanf(fp,"%d",&temp)!=EOF){
		DC_Shift+=temp;
		cnt++;
	}
	DC_Shift/=cnt;
	fclose(fp);

	//STEP2:NORMALISATION OF DATA
	
	//Now we have to process data so we open this file
	char filename[100],filename1[100];
	strcpy(filename,str);
	strcpy(filename1,"Normalised.txt");
	

	fp=fopen(filename,"r");
	//This is to store the normalised data
	fp1=fopen(filename1,"w");
	//For skipping the first 5 lines
	cnt=0;
	while(fgets(InitialStrings,80,fp)!=NULL && cnt<5){
		cnt++;
	}
    temp=0;
	//after that cool edit will add some headers at the beginning.we have to remove those
	cnt=0;
	while(fscanf(fp,"%d",&temp)!=EOF){
		if(temp<-1)
			break;
		if(temp>1)
			break;
		cnt++;
	}
	//printf("cnt %d",cnt);
	cnt=0;
	temp=0;
	long double max=INT_MIN;
	//position of file pointer
	long int position=ftell(fp);
	//To find max for normalisation
	while(fscanf(fp,"%d",&temp)!=EOF){
		if(temp<0){
			if(max<(DC_Shift-temp)){
				max=(DC_Shift-temp);
			}
		}
		else{
			if(max<(temp-DC_Shift)){
				max=temp-DC_Shift;
			}
		}
		cnt++;
	}
	//setting the file pointer position
	fseek(fp,position,SEEK_SET);
	temp=0;
	//storing the normalised data in file
	long double NormalisationFactor=((float)5000/max);
	//Normalising the Data and storing it in file
	while(fscanf(fp,"%d",&temp)!=EOF){
		long int t=(temp-DC_Shift)*NormalisationFactor;
		fprintf(fp1,"%d\n",t);
	}
	fclose(fp);
	fclose(fp1);

	/*STEP3: CALCULATION OF ENERGY AND ZCR FOR THE NORMALISED DATA BY USING
		     USING WINDOW SIZE OF 320  AND ALSO CALCULATING 
			 THE AMBIENT NOISE */

	
	//calculating Energy and zcr for the produced Normalised data and also AmbientNoise
	strcpy(filename,"EnergyZcr.txt");
	//strcat(filename,str);
	fp=fopen(filename1,"r");
	fp1=fopen(filename,"w");
	long int cnt1=0;
	position=0;
	while(1){
		cnt=0;
		temp=0;
		long int sum=0,zcr=0,pos=0,neg=0;
		long int buf[100]={0};
		while(fscanf(fp,"%d",&temp)!=EOF && cnt<100){
			buf[cnt]=(temp);
			sum+=buf[cnt]*buf[cnt];
			if(buf[cnt]>0)
				pos=1;
			if(buf[cnt]<0)
				neg=1;
			if(pos==1 && neg==1){
				zcr++;
				if(buf[cnt]>0)
					neg=0;
				if(buf[cnt]<0)
					pos=0;
			}
			cnt++;
		}
		if(fscanf(fp,"%d",&temp)==EOF){
			break;
		}
		long double energy=(long double)sum/100;
		fprintf(fp1,"%d ",cnt1);
		fprintf(fp1,"%g ",energy);
		fprintf(fp1,"%d ",zcr);
		fprintf(fp1,"%d\n",ftell(fp));
		//printf("energy:%g\n",energy);
		Ambient_noise+=energy;

		cnt1++;
	}
	//printf("count1:%d\n",cnt1);
	Ambient_noise/=cnt1;
	//printf("count1:%g\n",Ambient_noise);
	fclose(fp);
	fclose(fp1);
	
	
	
	/*STEP 4: FINDING THE 5 STEADY FRAMES*/

	
	
	fp=fopen(filename,"r");
	fp1=fopen(filename1,"r");
	cnt=0;
	//cnt1=0;
	temp=0;
	float tempEnergy=0;
	int tempZcr=0;
	int framePosition;
	int frameNo;
	int p;
	/*long double max=INT_MIN;
	int index;
	//printf("frame NOs\n");
	while(fscanf(fp,"%d %f %d %d",&temp,&tempEnergy,&tempZcr,&position)!=EOF){
		if(tempEnergy>max){
			max=tempEnergy;
			index=temp;
		}	
	}
	//printf("%g %d \n",max,index);
	//printf("%d\n",cnt1);
	*/
	fseek(fp,0,SEEK_SET);
	while(fscanf(fp,"%d %f %d %d",&temp,&tempEnergy,&tempZcr,&position)!=EOF){
			if(tempEnergy>(Ambient_noise)){
				if(cnt==0){
					framePosition=p;
					frameNo=temp;
					//printf("temp");
				}
				cnt++;
			}
			else{
				cnt=0;
			}
			if(cnt==5){
				break;
			}
			p=position;
	}
	fclose(fp);
	//printf("Frame no:%d\n",frameNo);
	strcpy(filename,"SteadyFrames.txt");
	//strcat(filename,str);
	fp=fopen(filename,"w");
	fseek(fp1,framePosition,SEEK_SET);
	cnt=1;
	while(fscanf(fp1,"%d",&temp)!=EOF && cnt<=85*100){
		fprintf(fp,"%d %d\n",cnt,temp);
		cnt++;
	}
	int No_of_frames_Needed=85-(cnt1-frameNo);
	//printf("frames Needed:%d\n",No_of_frames_Needed);
	while(No_of_frames_Needed>0){
		int temp1=cnt;
		fseek(fp1,0,SEEK_SET);
		//printf("cnt:%d\n",temp1);
		while(fscanf(fp1,"%d",&temp)!=EOF && cnt<=8500){
			//printf("%d\n",temp+(frameNo)*320);
			if(cnt==temp1+(frameNo)*100){
				//printf("***\n");
				break;
			}
			fprintf(fp,"%d %d\n",cnt,temp);
			cnt++;
		}
		No_of_frames_Needed=No_of_frames_Needed-frameNo;
	}
	//printf("%d\n",cnt);
	fclose(fp);
	fclose(fp1);
    

}


//This Function is used to find the cepstral coefficients
void CalculateRAC(){
	FILE *fp;
	//Opening test.txt in read mode
	char filename[100];
	strcpy(filename,"SteadyFrames.txt");
	fp=fopen(filename,"r");
	//fp1=fopen("testCi.txt","w");
	int count=0;
	//global_i=0;
	while(count<85){
		count++;
		double test[100]={0};
		int temp=0,t=0;
		int cnt=0;
		//Reading Frames from the SteadyFrames.txt file.
		while(cnt<100 && fscanf(fp,"%d %d",&t,&temp)!=EOF){
			test[cnt++]=temp;
		}
		temp=0;
		//calculation of Ri's
		double R[13];
		for(int i=0;i<=12;i++){
			R[i]=0;
			for(int j=0;j<100-i;j++){
				R[i]=R[i]+(test[j]*test[j+i]);
			}
		}
		
		//Calculation of Ai's using durbins algorithm
		double E[13];
		double a[13][13];
		double k[13];
		double A[12];
		double C[13];
		E[0]=R[0];
		for(int i=1;i<=12;i++){
			double sum=0;
			for(int j=1;j<=i-1;j++){
				sum=sum+(a[i-1][j]*R[i-j]);
			}
			k[i]=(R[i]-sum)/E[i-1];
			a[i][i]=k[i];
			for(int j=1;j<=i-1;j++){
				a[i][j]=a[i-1][j]-(k[i]*a[i-1][i-j]);
			}
			E[i]=(1-(k[i]*k[i]))*E[i-1];
		}
		for(int i=1;i<=12;i++){
			A[i-1]=a[12][i];
		}
		
		//Calculation of Ci's
		double w=0;
		C[0]=log(R[0]*R[0]);
		for(int i=1;i<=12;i++){
			double sum=0;
			double m=i;
			for(int j=1;j<=i-1;j++){
				double k=j; 
				sum=sum+((k*C[j]*A[i-j-1])/m);
			}
			//printf("w:%lf ",w);
			C[i]=(A[i-1]+sum);
			//printf("C[i]:%lf\n",C[global_i][i]);
		}
		for(int i=1;i<=12;i++){
			double q=12;
			C[i]=C[i]*(1+(q/2)*sin((pi*i)/q));
			
		}
		//fprintf(fp1,"\n");
		for(int i=1;i<=12;i++){
			Universe[global_i][i]=C[i];
		}
		global_i++;
		
	}
	fclose(fp);
	//fclose(fp1);
}

















//This Function generates the observation sequence
void Observation_Sequence_Generation(){

	long double tokhurasWeights[12]={1.0,3.0,7.0,13.0,19.0,22.0,25.0,33.0,42.0,50.0,56.0,61.0};
		
	//Calculating Tokhuras Distance and Predicting the Output
	for(int i=1;i<=85;i++){
		long double min=INT_MAX;
		int min_index=0;
		for(int j=1;j<=32;j++){
			long double sum=0;
			for(int k=1;k<=12;k++){
				long double diff=Universe[i][k]-codebook[j][k];
				sum=sum+tokhurasWeights[k-1]*diff*diff;
			}
			if(sum<min){
				min=sum;
				min_index=j;
			}
		}
		O[i]=min_index;
		
	}
}













//Feed Forward Model or Bakis Model
void Feed_Forward(){
	//Initializing Pi values
	for(int i=1;i<=N;i++){
		if(i==1)
			Pi[i]=1.0;
		else
			Pi[i]=0.0;
	}

	//Initializing A values
	for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++){
			if(i==j)
				A[i][j]=0.8;
			else if(i+1==j)
				A[i][j]=0.2;
			else
				A[i][j]=0.0;
		}
	}
	
	//Initializing B values
	for(int i=1;i<=N;i++){
		for(int j=1;j<=M;j++){
			B[i][j]=0.03125;
		}
	}
}


//Forward Procedure
void Forward_Procedure(){
	//Initialization
	for(int i=1;i<=N;i++){
		Alpha[1][i]=Pi[i]*B[i][O[1]];
		
	}

	//Induction
	for(int t=2;t<=T;t++){
		for(int j=1;j<=N;j++){
			long double s=0;			
			for(int i=1;i<=N;i++){
				s=s+(Alpha[t-1][i]*A[i][j]);
			}
			Alpha[t][j]=s*B[j][O[t]];		
		}
	}
	
	//Termination
	long double Prob=0;
	for(int i=1;i<=N;i++){
		Prob=Prob+Alpha[T][i];
		
	}
	
}


//Backward_procedure
void Backward_Procedure(){
	//Initialization
	for(int i=1;i<=N;i++){
		Beta[T][i]=1.0;
	}

	//Induction
	for(int t=T-1;t>0;t--){
		for(int i=1;i<=N;i++){
			long double s=0;
			for(int j=1;j<=N;j++){
				s=s+(A[i][j]*B[j][O[t+1]]*Beta[t+1][j]);
			}
			Beta[t][i]=s;
		}
	}
}



//Viterbi Test
long double Viterbi_test(int p){
	//Initialization
	for(int i=1;i<=N;i++){
		Delta[1][i]=AVG_Pi[p][i]*AVG_B[p][i][O[1]];
		shi[1][i]=0;
	}

	//Recursion
	for(int t=2;t<=T;t++){
		for(int j=1;j<=N;j++){
			long double temp=0;
			int index=0;
			for(int i=1;i<=N;i++){
				if(i==1){
					temp=Delta[t-1][i]*AVG_A[p][i][j];
					index=i;
				}
				else if(temp<(Delta[t-1][i]*AVG_A[p][i][j])){
					temp=Delta[t-1][i]*AVG_A[p][i][j];
					index=i;					
				}
			}
			Delta[t][j]=temp*AVG_B[p][j][O[t]];
			shi[t][j]=index;
		}
	}
	//Termination
	long double temp=0;
	long double index=0;

	for(int i=1;i<=N;i++){
		if(i==1){
			temp=Delta[T][i];
			index=i;
		}
		else if(temp<(Delta[T][i])){
			temp=Delta[T][i];
			index=i;					
		}
	}
	long double P_Star=temp;
	int Q_Star[T+1];
	Q_Star[T]=index;

	//Back Tracking
	for(int t=T-1;t>0;t--){
		Q_Star[t]=shi[t+1][Q_Star[t+1]];
	}
	return P_Star;
}





//Viterbi
void Viterbi(){
	//Initialization
	for(int i=1;i<=N;i++){
		Delta[1][i]=Pi[i]*B[i][O[1]];
		shi[1][i]=0;
	}

	//Recursion
	for(int t=2;t<=T;t++){
		for(int j=1;j<=N;j++){
			long double temp=0;
			int index=0;
			for(int i=1;i<=N;i++){
				if(i==1){
					temp=Delta[t-1][i]*A[i][j];
					index=i;
				}
				else if(temp<(Delta[t-1][i]*A[i][j])){
					temp=Delta[t-1][i]*A[i][j];
					index=i;					
				}
			}
			Delta[t][j]=temp*B[j][O[t]];
			shi[t][j]=index;
		}
	}
	//Termination
	long double temp=0;
	long double index=0;

	for(int i=1;i<=N;i++){
		if(i==1){
			temp=Delta[T][i];
			index=i;
		}
		else if(temp<(Delta[T][i])){
			temp=Delta[T][i];
			index=i;					
		}
	}
	long double P_Star=temp;
	int Q_Star[T+1];
	Q_Star[T]=index;

	//Back Tracking
	for(int t=T-1;t>0;t--){
		Q_Star[t]=shi[t+1][Q_Star[t+1]];
	}
	//Writing Pstar to the file
	fprintf(globalFile,"P_star:%g\n",P_Star);
	fprintf(globalFile,"Q_Star Values: \n");
	for(int t=1;t<=T;t++){
		fprintf(globalFile,"%d ",Q_Star[t]);
	}
	fprintf(globalFile,"\n");
	
}


//Gamma Procedure
void Gamma_Procedure(){
	
	//Calculating Gamma Values using Alpha and Beta
	for(int t=1;t<=T;t++){
		long double s=0;
		for(int i=1;i<=N;i++){
			s=s+Alpha[t][i]*Beta[t][i];		
		}
		for(int i=1;i<=N;i++){
			Gamma[t][i]=(Alpha[t][i]*Beta[t][i])/s;
		}
	}

}


//Baum_Welch Function
void Baum_Welch(){

	//Finding Eta values
	for(int t=1;t<T;t++){
		long double s=0;
		for(int i=1;i<=N;i++){
			for(int j=1;j<=N;j++){
				s=s+(Alpha[t][i]*A[i][j]*B[j][O[t+1]]*Beta[t+1][j]);
			}
		}
		for(int i=1;i<=N;i++){
			for(int j=1;j<=N;j++){
				Eta[t][i][j]=(Alpha[t][i]*A[i][j]*B[j][O[t+1]]*Beta[t+1][j])/s;
			}
		}
	}

}



//Re_Estimation or Updation of Model
void Re_Estimation(){

	//Updating Pi values
	for(int i=1;i<=N;i++){
		Pi[i]=Gamma[1][i];
	}
	
	//Updating A Values
	for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++){
			long double s=0;
			for(int t=1;t<T;t++){
					s=s+Gamma[t][i];
			}
			long double s1=0;
			for(int t=1;t<T;t++){
				s1=s1+Eta[t][i][j];
			}
			A[i][j]=s1/s;
		}
	}

	//Updating B values
	for(int i=1;i<=N;i++){
		for(int j=1;j<=M;j++){
			long double s=0;
			for(int t=1;t<=T;t++){
					s=s+Gamma[t][i];
			}
			long double s1=0;
			for(int t=1;t<=T;t++){
				if(O[t]==j){
					s1=s1+Gamma[t][i];
				}
			}
			B[i][j]=s1/s;
			if(B[i][j]==0){
				B[i][j]=1.0e-30;
			}
		}
	}
	

}

//Stochastic Check function to maintain the row sum to be 1 in A and B matrices.
void Stochastic_Check(){

	//For Maintaing A matrix row sum to be 1
	for(int i=1;i<=N;i++){
		long double s=0;
		long double max=INT_MIN;
		int index=0;
		for(int j=1;j<=N;j++){
			s=s+A[i][j];
			if(max<A[i][j]){
				max=A[i][j];
				index=j;
			}
		}
		if(s!=1){
			if(s>1){
				A[i][index]=A[i][index]-(s-1);
			}
			else if(s<1){
				A[i][index]=A[i][index]+(1-s);
			}
		}
	}

	//Maintaing B Matrix row sum to be 1
	for(int i=1;i<=N;i++){
		long double s=0;
		long double max=INT_MIN;
		int index=0;
		for(int j=1;j<=M;j++){
			s=s+B[i][j];
			if(max<B[i][j]){
				max=B[i][j];
				index=j;
			}
		}
		if(s!=1){
			if(s>1){
				B[i][index]=B[i][index]-(s-1);
			}
			else if(s<1){
				B[i][index]=B[i][index]+(1-s);
			}
		}
	}
}


//Stochastic check on the Averaged A and B values
void Stochastic_Check_AVG(int p){
	//For maintaing averaged A values row sum to be 1
	for(int i=1;i<=N;i++){
		long double s=0;
		long double max=INT_MIN;
		int index=0;
		for(int j=1;j<=N;j++){
			s=s+AVG_A[p][i][j];
			if(max<AVG_A[p][i][j]){
				max=AVG_A[p][i][j];
				index=j;
			}
		}
		if(s!=1){
			if(s>1){
				AVG_A[p][i][index]=AVG_A[p][i][index]-(s-1);
			}
			else if(s<1){
				AVG_A[p][i][index]=AVG_A[p][i][index]+(1-s);
			}
		}
	}

	//For maintaining Averaged B values row sum to be 1
	for(int i=1;i<=N;i++){
		long double s=0;
		long double max=INT_MIN;
		int index=0;
		for(int j=1;j<=M;j++){
			s=s+AVG_B[p][i][j];
			if(max<AVG_B[p][i][j]){
				max=AVG_B[p][i][j];
				index=j;
			}
		}
		if(s!=1){
			if(s>1){
				AVG_B[p][i][index]=AVG_B[p][i][index]-(s-1);
			}
			else if(s<1){
				AVG_B[p][i][index]=AVG_B[p][i][index]+(1-s);
			}
		}
	}
}



void Solution_To_Problem1(){
	Forward_Procedure();
	Backward_Procedure();
}


void Solution_To_Problem2(){
	Viterbi();
}


void Solution_To_Problem3(){
	Baum_Welch();
	Gamma_Procedure();
	Re_Estimation();
}

void train(){
	char fileA[80];
	strcpy(fileA,Names[index-1]);
	char fileB[80];
	strcpy(fileB,Names[index-1]);
	char filePi[80];
	strcpy(filePi,Names[index-1]);
	strcat(fileA,"A.txt");
	strcat(fileB,"B.txt");
	strcat(filePi,"Pi.txt");
	FILE *fp=fopen(filePi,"w");
	FILE *fp1=fopen(fileA,"w");
	FILE *fp2=fopen(fileB,"w");
	FILE *fp3=fopen("codebook.txt","w");
	LBG();
	//cout<<global_i<<endl;
	//Printing the CodeBook
	for(int i=1;i<=32;i++){
		for(int j=1;j<=12;j++){
			printf("%g ",codebook[i][j]);
			fprintf(fp3,"%g ",codebook[i][j]);
		}
		printf("\n");
		fprintf(fp3,"\n");
	}
	fclose(fp3);
	//Training the 20 recorded files
	int cl=0;
	while(cl<3){
		for(int i=index-1;i<index;i++){
			for(int l=1;l<=N;l++){
				Pi[l]=AVG_Pi[i][l];
				AVG_Pi[i][l]=0;
			}
			for(int l=1;l<=N;l++){
				for(int m=1;m<=N;m++){
					A[l][m]=AVG_A[i][l][m];
					AVG_A[i][l][m]=0;
				
				}
				
			}
			
			for(int l=1;l<=N;l++){
				for(int m=1;m<=M;m++){
					B[l][m]=AVG_B[i][l][m];
					AVG_B[i][l][m]=0;
				
				}
				
			}
			for(int j=0;j<20;j++){
				global_i=1;
				char s[100];
				strcpy(s,Names[i]);
				int len=strlen(s);
				if(j<10){
					s[len]='0'+j;
					s[len+1]='\0';			
				}
				else if(j<20){
					s[len]='1';
					s[len+1]='0'+(j%10);
					s[len+2]='\0';				
				}
				strcat(s,".txt");
				//richTextBox1->Text+=gcnew String(s);
				//richTextBox1->Text+=gcnew String("\n");	
				//Debug.WriteLine(s);
				//Debug.WriteLine("\n");
				CalculateSteadyFrames(s);
				CalculateRAC();
				//LBG();
				Observation_Sequence_Generation();
				Feed_Forward();
				int count=0;
				while(count<200){
					Solution_To_Problem1();
					Solution_To_Problem2();
					Solution_To_Problem3();
					Stochastic_Check();
					count++;
				}
				//printf("Pi values %d %d:\n",i,j);
				for(int l=1;l<=N;l++){
					AVG_Pi[i][l]+=Pi[l];
					//printf("%g ",Pi[l]);
				}
				//printf("\n");
				//printf("A values:\n");
				for(int l=1;l<=N;l++){
					for(int m=1;m<=N;m++){
						AVG_A[i][l][m]+=A[l][m];
						//printf("%g ",A[l][m]);
					}
					//printf("\n");
				}
				//printf("\n");
				//printf("B values:\n");
				for(int l=1;l<=N;l++){
					for(int m=1;m<=M;m++){
						AVG_B[i][l][m]+=B[l][m];
						//printf("%g ",B[l][m]);
					}
					//printf("\n");
				}
				//printf("\n");
			
			}
			//printf("Average Values For %d \n:",i);
			//printf("Final Pi values:\n");
			for(int l=1;l<=N;l++){
				AVG_Pi[i][l]/=20;
				//printf("%g ",AVG_Pi[i][l]);
			}
			//printf("\nFinal A values:\n");
			for(int l=1;l<=N;l++){
				for(int m=1;m<=N;m++){
					AVG_A[i][l][m]/=20;
					//printf("%g ",AVG_A[i][l][m]);
				}
				//printf("\n");
			}
			//printf("\nFinal B values:\n");
			for(int l=1;l<=N;l++){
				for(int m=1;m<=M;m++){
					AVG_B[i][l][m]/=20;
					//printf("%g ",AVG_B[i][l][m]);
				}
				//printf("\n");
			}
			Stochastic_Check_AVG(i);
			fprintf(globalFile,"Pi Values:\n");
			for(int l=1;l<=N;l++){
				fprintf(globalFile,"%g ",AVG_Pi[i][l]);
			}
			fprintf(globalFile,"\nA Values:\n");
			for(int l=1;l<=N;l++){
				for(int m=1;m<=N;m++){
					//AVG_A[i][l][m]+=A[l][m];
					fprintf(globalFile,"%g ",AVG_A[i][l][m]);
				}
				fprintf(globalFile,"\n");
			}
			fprintf(globalFile,"\nB Values:\n");
			for(int l=1;l<=N;l++){
				for(int m=1;m<=M;m++){
					//AVG_A[i][l][m]+=A[l][m];
					fprintf(globalFile,"%g ",AVG_B[i][l][m]);
				}
				fprintf(globalFile,"\n");
			}			
			fprintf(globalFile,"------------------------------------------------------------------------------------\n");
					
		}
		cl++;
	}
	for(int i=index-1;i<index;i++){
		printf("Average Values For %d: \n",i);
		fprintf(globalFile,"Average Values For %d: \n",i);
		printf("Final Pi values:\n");
		fprintf(globalFile,"Final Pi values:\n");
		for(int l=1;l<=N;l++){
			//AVG_Pi[i][l]/=20;
			printf("%g ",AVG_Pi[i][l]);
			fprintf(globalFile,"%g ",AVG_Pi[i][l]);
			fprintf(fp,"%g ",AVG_Pi[i][l]);
		}
		printf("\nFinal A values:\n");
		fprintf(fp,"\n");
		fprintf(globalFile,"\nFinal A values:\n");
		for(int l=1;l<=N;l++){
			for(int m=1;m<=N;m++){
				//AVG_A[i][l][m]/=20;
				printf("%g ",AVG_A[i][l][m]);
				fprintf(globalFile,"%g ",AVG_A[i][l][m]);
				fprintf(fp1,"%g ",AVG_A[i][l][m]);
			}
			printf("\n");
			fprintf(globalFile,"\n");
			fprintf(fp1,"\n");
		}
		fprintf(fp1,"\n");
		printf("\nFinal B values:\n");
		fprintf(globalFile,"\nFinal B values:\n");
		for(int l=1;l<=N;l++){
			for(int m=1;m<=M;m++){
				//AVG_B[i][l][m]/=20;
				printf("%g ",AVG_B[i][l][m]);
				fprintf(globalFile,"%g ",AVG_B[i][l][m]);
				fprintf(fp2,"%g ",AVG_B[i][l][m]);
			}
			printf("\n");
			fprintf(globalFile,"\n");
			fprintf(fp2,"\n");
		}
		fprintf(fp2,"\n");
		printf("--------------------------------------------------------------------\n");
		fprintf(globalFile,"------------------------------------------------------------------------------------\n");
	}
	fclose(fp);
	fclose(fp1);
	fclose(fp2);

}





int test(){
	//int ch=1;
	//while(ch==1){
		system("Recording_Module.exe 3 input.wav input.txt");
		global_i=1;
		CalculateSteadyFrames("input.txt");
		CalculateRAC();
		//LBG();
		Observation_Sequence_Generation();
		long double max=INT_MIN;
		int predicted_output=-1;
		FILE *fp=fopen("testing.txt","w");
		for(int i=0;i<index;i++){
			for(int l=1;l<=N;l++){
				fprintf(fp,"%g ",AVG_Pi[i][l]);
			}
			fprintf(fp,"\n");
			for(int l=1;l<=N;l++){
				for(int m=1;m<=N;m++){
					fprintf(fp,"%g ",AVG_A[i][l][m]);
				}
				fprintf(fp,"\n");
			}
			fprintf(fp,"\n");
			for(int l=1;l<=N;l++){
				for(int m=1;m<=M;m++){
					fprintf(fp,"%g ",AVG_B[i][l][m]);
					
				}
				fprintf(fp,"\n");
			}
			fprintf(fp,"\n");

		}
		for(int k=0;k<index;k++){
			long double prob=Viterbi_test(k);
			printf("%d %g \n",k,prob);
			if(prob>max){
				max=prob;
				predicted_output=k;
			}
		}
		printf("The Pridicted Output for this File is %s. \n",Names[predicted_output]);
		//cout<<"Enter 1 for testing agiain and 0 not to test:"<<endl;
	return predicted_output;
}

namespace Voiced_Based_ContactSearch {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Summary for Form1
	/// </summary>
	public ref class Form1 : public System::Windows::Forms::Form
	{
	public:
		Form1(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~Form1()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Label^  label1;
	protected: 
	private: System::Windows::Forms::TextBox^  textBox1;
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::TextBox^  textBox2;
	private: System::Windows::Forms::Button^  button1;

	private: System::Windows::Forms::Button^  button3;
	private: System::Windows::Forms::Label^  label3;
	private: System::Windows::Forms::RichTextBox^  richTextBox1;




	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(Form1::typeid));
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->textBox2 = (gcnew System::Windows::Forms::TextBox());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->button3 = (gcnew System::Windows::Forms::Button());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->richTextBox1 = (gcnew System::Windows::Forms::RichTextBox());
			this->SuspendLayout();
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->BackColor = System::Drawing::Color::Cornsilk;
			this->label1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label1->Location = System::Drawing::Point(51, 57);
			this->label1->Margin = System::Windows::Forms::Padding(4, 0, 4, 0);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(60, 20);
			this->label1->TabIndex = 0;
			this->label1->Text = L"Name:";
			this->label1->Click += gcnew System::EventHandler(this, &Form1::label1_Click);
			// 
			// textBox1
			// 
			this->textBox1->BackColor = System::Drawing::SystemColors::Window;
			this->textBox1->Location = System::Drawing::Point(163, 57);
			this->textBox1->Margin = System::Windows::Forms::Padding(4, 5, 4, 5);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(290, 26);
			this->textBox1->TabIndex = 1;
			this->textBox1->TextChanged += gcnew System::EventHandler(this, &Form1::textBox1_TextChanged);
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->BackColor = System::Drawing::Color::Cornsilk;
			this->label2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->label2->Location = System::Drawing::Point(51, 109);
			this->label2->Margin = System::Windows::Forms::Padding(4, 0, 4, 0);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(99, 20);
			this->label2->TabIndex = 2;
			this->label2->Text = L"ContactNo:";
			this->label2->Click += gcnew System::EventHandler(this, &Form1::label2_Click);
			// 
			// textBox2
			// 
			this->textBox2->Location = System::Drawing::Point(163, 106);
			this->textBox2->Margin = System::Windows::Forms::Padding(4, 5, 4, 5);
			this->textBox2->Name = L"textBox2";
			this->textBox2->Size = System::Drawing::Size(290, 26);
			this->textBox2->TabIndex = 3;
			this->textBox2->TextChanged += gcnew System::EventHandler(this, &Form1::textBox2_TextChanged);
			// 
			// button1
			// 
			this->button1->BackColor = System::Drawing::SystemColors::Desktop;
			this->button1->ForeColor = System::Drawing::SystemColors::ButtonFace;
			this->button1->Location = System::Drawing::Point(234, 154);
			this->button1->Margin = System::Windows::Forms::Padding(4, 5, 4, 5);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(112, 47);
			this->button1->TabIndex = 4;
			this->button1->Text = L"Insert";
			this->button1->UseVisualStyleBackColor = false;
			this->button1->Click += gcnew System::EventHandler(this, &Form1::button1_Click);
			// 
			// button3
			// 
			this->button3->BackColor = System::Drawing::SystemColors::Desktop;
			this->button3->FlatStyle = System::Windows::Forms::FlatStyle::Flat;
			this->button3->ForeColor = System::Drawing::SystemColors::ButtonFace;
			this->button3->Location = System::Drawing::Point(671, 92);
			this->button3->Margin = System::Windows::Forms::Padding(4, 5, 4, 5);
			this->button3->Name = L"button3";
			this->button3->Size = System::Drawing::Size(112, 43);
			this->button3->TabIndex = 6;
			this->button3->Text = L"Search";
			this->button3->UseVisualStyleBackColor = false;
			this->button3->Click += gcnew System::EventHandler(this, &Form1::button3_Click);
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 15.75F, static_cast<System::Drawing::FontStyle>((System::Drawing::FontStyle::Bold | System::Drawing::FontStyle::Italic)), 
				System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(0)));
			this->label3->Location = System::Drawing::Point(52, 14);
			this->label3->Margin = System::Windows::Forms::Padding(4, 0, 4, 0);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(313, 25);
			this->label3->TabIndex = 7;
			this->label3->Text = L"Voice Based Contact Search";
			// 
			// richTextBox1
			// 
			this->richTextBox1->Location = System::Drawing::Point(18, 238);
			this->richTextBox1->Margin = System::Windows::Forms::Padding(4, 5, 4, 5);
			this->richTextBox1->Name = L"richTextBox1";
			this->richTextBox1->Size = System::Drawing::Size(846, 176);
			this->richTextBox1->TabIndex = 8;
			this->richTextBox1->Text = L"";
			this->richTextBox1->TextChanged += gcnew System::EventHandler(this, &Form1::richTextBox1_TextChanged);
			// 
			// Form1
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(9, 20);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->BackColor = System::Drawing::Color::Cornsilk;
			this->BackgroundImage = (cli::safe_cast<System::Drawing::Image^  >(resources->GetObject(L"$this.BackgroundImage")));
			this->BackgroundImageLayout = System::Windows::Forms::ImageLayout::Stretch;
			this->ClientSize = System::Drawing::Size(884, 435);
			this->Controls->Add(this->richTextBox1);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->button3);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->textBox2);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->textBox1);
			this->Controls->Add(this->label1);
			this->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point, 
				static_cast<System::Byte>(0)));
			this->Margin = System::Windows::Forms::Padding(4, 5, 4, 5);
			this->Name = L"Form1";
			this->Text = L"Form1";
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void label1_Click(System::Object^  sender, System::EventArgs^  e) {
			 }
private: System::Void label2_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void textBox1_TextChanged(System::Object^  sender, System::EventArgs^  e) {
		
}
private: System::Void textBox2_TextChanged(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {	
		initial();
		string s = marshal_as<string>(textBox1->Text);
		FILE *fp,*fp1,*fp2;
		fp=fopen("Names.txt","a");
		fp1=fopen("ContactNo.txt","a");
		if (s.length() < sizeof(Names[index])) 
			sprintf(Names[index], "%s", s.c_str());
		s = marshal_as<string>(textBox2->Text);
		if (s.length() < sizeof(ContactNo[index])) 
			sprintf(ContactNo[index], "%s", s.c_str());
		fprintf(fp,"%s\n",Names[index]);
		fprintf(fp1,"%s\n",ContactNo[index]);
		fclose(fp);
		fclose(fp1);
		index++;
		int count=0;
		cout<<"Record the Name for 20 times(1 time per Recording)"<<endl;
		while(count<20){
			char temp[1000]="Recording_Module.exe 3 ";
			char temp1[100]="";
			strcat(temp,Names[index-1]);
			strcat(temp1,Names[index-1]);
			if(count<10){
				int l=strlen(temp);
				int l1=strlen(temp1);
				temp[l]='0'+count;
				temp[l+1]='\0';
				temp1[l1]='0'+count;
				temp1[l1+1]='\0';
			}
			else if(count<20){
				int l=strlen(temp);
				temp[l]='1';
				temp[l+1]='0'+(count%10);
				temp[l+2]='\0';
				int l1=strlen(temp1);
				temp1[l1]='1';
				temp1[l1+1]='0'+(count%10);
				temp1[l1+2]='\0';
			}
			strcat(temp,".wav ");
			strcat(temp,Names[index-1]);
			if(count<10){
				int l=strlen(temp);
				temp[l]='0'+count;
				temp[l+1]='\0';
			}
			else if(count<20){
				int l=strlen(temp);
				temp[l]='1';
				temp[l+1]='0'+(count%10);
				temp[l+2]='\0';
			}
			strcat(temp,".txt");
			strcat(temp1,".wav");
			richTextBox1->Text+=gcnew String(temp);
			richTextBox1->Text+=gcnew String("\n");
			richTextBox1->Text+=gcnew String(temp1);
			richTextBox1->Text+=gcnew String("\n");
			system(temp);	
			count++;
		}	
		richTextBox1->Text+=gcnew String("Please Wait Until Traning gets completed message \n");
		global_i=1;
		Gflag=0;
		for(int i=0;i<index;i++){
			for(int j=0;j<20;j++){
				char s[100];
				strcpy(s,Names[i]);
				int len=strlen(s);
				if(j<10){
				s[len]='0'+j;
				s[len+1]='\0';			
				}
				else if(j<20){
					s[len]='1';
					s[len+1]='0'+(j%10);
					s[len+2]='\0';				
				}
				strcat(s,".txt");
				CalculateSteadyFrames(s);
				CalculateRAC();
			}
		}
		global_i--;
		richTextBox1->Text+=gcnew String("Traning Started:\n");
		for(int i=0;i<index-1;i++){
			char fileA[80];
			strcpy(fileA,Names[i]);
			char fileB[80];
			strcpy(fileB,Names[i]);
			char filePi[80];
			strcpy(filePi,Names[i]);
			strcat(fileA,"A.txt");
			strcat(fileB,"B.txt");
			strcat(filePi,"Pi.txt");
			fp=fopen(filePi,"r");
			fp1=fopen(fileA,"r");
			fp2=fopen(fileB,"r");
			for(int l=1;l<=N;l++){
				fscanf(fp,"%Lf",&AVG_Pi[i][l]);
			}
			for(int l=1;l<=N;l++){
				for(int m=1;m<=N;m++){
					fscanf(fp1,"%Lf",&AVG_A[i][l][m]);
				}
			}
			for(int l=1;l<=N;l++){
				for(int m=1;m<=M;m++){
					fscanf(fp2,"%Lf",&AVG_B[i][l][m]);
					
				}
			}
			fclose(fp);
			fclose(fp1);
			fclose(fp2);
		}
		train();
		richTextBox1->Text=gcnew String("Insertion Completed\n");
		textBox1->Text=gcnew String("");
		textBox2->Text=gcnew String("");
		wstring ws=s2ws("insertion_successfull.wav");
		LPCWSTR r=ws.c_str();
		PlaySound(r, NULL, SND_FILENAME);

}
/*private: System::Void button2_Click(System::Object^  sender, System::EventArgs^  e) {
	if(index==0){
		richTextBox1->Text=gcnew String("Deleted Successfully\n");
		textBox1->Text=gcnew String("");
		textBox2->Text=gcnew String("");
	}
	else{
		char temp[1000],temp1[1000];
		bool flag=0;
		FILE *fp=fopen("Deleted.txt","a");
		string s = marshal_as<string>(textBox1->Text);
		if (s.length() < sizeof(temp)) 
			sprintf(temp, "%s", s.c_str());
		s = marshal_as<string>(textBox2->Text);
		if (s.length() < sizeof(temp1)) 
			sprintf(temp1, "%s", s.c_str());
		fprintf(fp,"%s\n",temp);
		for(int i=0;i<index;i++){
			if(strcmp(temp,Names[i])==0 && strcmp(temp1,ContactNo[index])==0){
				flag=1;
			}
			if(flag==1 && i!=index-1){
				strcpy(Names[i],Names[i+1]);
				strcpy(ContactNo[i],ContactNo[i+1]);
			}
		}
		index=index-1;
		Gflag=1;
		richTextBox1->Text=gcnew String("Deleted Successfully\n");
		textBox1->Text=gcnew String("");
		textBox2->Text=gcnew String("");
	}		 
}*/
private: System::Void button3_Click(System::Object^  sender, System::EventArgs^  e) {
		if(Gflag==1){
			initial();
			FILE *fp,*fp1,*fp2;
			//char c='0'+index;
			//richTextBox1->Text+=gcnew String(c);
			//fp3=fopen("testing.txt","w");
			for(int i=0;i<index;i++){
				//richTextBox1->Text+=gcnew String(Names[i]);
				//richTextBox1->Text+=gcnew String("\n");
				char fileA[80];
				strcpy(fileA,Names[i]);
				char fileB[80];
				strcpy(fileB,Names[i]);
				char filePi[80];
				strcpy(filePi,Names[i]);
				strcat(fileA,"A.txt");
				//richTextBox1->Text+=gcnew String(fileA);
				//richTextBox1->Text+=gcnew String("\n");
				strcat(fileB,"B.txt");
				//richTextBox1->Text+=gcnew String(fileB);
				//richTextBox1->Text+=gcnew String("\n");
				strcat(filePi,"Pi.txt");
				//richTextBox1->Text+=gcnew String(filePi);
				//richTextBox1->Text+=gcnew String("\n");
				fp=fopen(filePi,"r");
				fp1=fopen(fileA,"r");
				fp2=fopen(fileB,"r");
				for(int l=1;l<=N;l++){
					fscanf(fp,"%Lf",&AVG_Pi[i][l]);
					//fprintf(fp3,"%Lf ",AVG_Pi[i][l]);
				}
				//fprintf(fp3,"\n");
				for(int l=1;l<=N;l++){
					for(int m=1;m<=N;m++){
						fscanf(fp1,"%Lf",&AVG_A[i][l][m]);
						//fprintf(fp3,"%Lf ",AVG_A[i][l][m]);
					}
					//fprintf(fp3,"\n");
				}
				//fprintf(fp3,"\n");
				for(int l=1;l<=N;l++){
					for(int m=1;m<=M;m++){
						fscanf(fp2,"%Lf",&AVG_B[i][l][m]);
						//fprintf(fp3,"%g ",AVG_B[i][l][m]);
					}
					//fprintf(fp3,"\n");
				}
				//fprintf(fp3,"\n");
				fclose(fp);
				fclose(fp1);
				fclose(fp2);
			}
		}
		Gflag=0;
		int op=test();
		richTextBox1->Text=gcnew String("Here Are Our Search Results:\n");
		richTextBox1->Text+=gcnew String("Name:");
		richTextBox1->Text+=gcnew String(Names[op]);
		richTextBox1->Text+=gcnew String("\n");
		richTextBox1->Text+=gcnew String("Contact Number:");
		richTextBox1->Text+=gcnew String(ContactNo[op]);
		richTextBox1->Text+=gcnew String("\n");
		wstring ws=s2ws("here_are_our_search_results.wav");
		LPCWSTR r=ws.c_str();
		PlaySound(r, NULL, SND_FILENAME);
		ws=s2ws("name_of_the_contact.wav");
		r=ws.c_str();
		PlaySound(r, NULL, SND_FILENAME);

		char temp1[100];
		strcpy(temp1,Names[op]);
		strcat(temp1,"0.wav");
		ws=s2ws(temp1);
		r=ws.c_str();
		//richTextBox1->Text+=gcnew String(temp1);
		PlaySound(r, NULL, SND_FILENAME);
		ws=s2ws("contact_number.wav");
		r=ws.c_str();
		PlaySound(r, NULL, SND_FILENAME);
		for(int i=0;i<strlen(ContactNo[op]);i++){
			temp1[0]=ContactNo[op][i];
			temp1[1]='\0';
			strcat(temp1,".wav");
			ws=s2ws(temp1);
			r=ws.c_str();
			//richTextBox1->Text+=gcnew String(temp1);
			PlaySound(r, NULL, SND_FILENAME);
		}
}
private: System::Void richTextBox1_TextChanged(System::Object^  sender, System::EventArgs^  e) {
}
private: System::Void label4_Click(System::Object^  sender, System::EventArgs^  e) {
		 }
private: System::Void button4_Click(System::Object^  sender, System::EventArgs^  e) {
		 
}
private: System::Void textBox3_TextChanged(System::Object^  sender, System::EventArgs^  e) {
		 }
};
}
//fclose(globalFile);
