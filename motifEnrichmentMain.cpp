#include <iostream>
#include <string>
#include <fstream>
//#include <stdlib>
#include "StringUtil.h"
#include "Nucleic.h";
#include <limits.h>
#include <pty.h>  /* for openpty and forkpty */
#include "py_hypergeom.h"
#include <stdio.h>
#include <map>
#include <vector>
#include <math.h>
#include <set>
#include <algorithm>
#define USE_JAVA_HYPERGEOMETRIC /*undefine this to use c++ hypergeometric function*/
using namespace std;

#ifdef USE_JAVA_HYPERGEOMETRIC
#define _BUFFER_SIZE 1024*1024*50
#define LOGMUL 1000
#define DEFEAULT_PSEUDO_COUNT 1
#define DEFAULT_N_PENALTY -1000000
#define A_J 0
#define C_J 1
#define T_J 2
#define G_J 3
#define NT_ALPHABET_SIZE 4
#define ALPHABET_STRING "ACTG"

class PSSM{
public:
	int **M;
	int length;
	int pseudoCount;
	int backgroundlogk[NT_ALPHABET_SIZE];
	int N_penalty;
	int numSeq;

	void setBackground(double fA,double fC,double fT,double fG)
	{
		double sumF=fA+fC+fT+fG+pseudoCount*NT_ALPHABET_SIZE;
		backgroundlogk[A_J]=int(LOGMUL*log((fA+pseudoCount)/sumF));
		backgroundlogk[C_J]=int(LOGMUL*log((fC+pseudoCount)/sumF));
		backgroundlogk[T_J]=int(LOGMUL*log((fT+pseudoCount)/sumF));
		backgroundlogk[G_J]=int(LOGMUL*log((fG+pseudoCount)/sumF));
	}

	void setBackground(const string& seq)
	{
		for(int i=0;i<NT_ALPHABET_SIZE;i++)
			backgroundlogk[i]=0; //this is for storing frequency

		int lseq=seq.length();
		for(int i=0;i<lseq;i++)
		{
			switch(seq[i])
			{
				case 'A':
					backgroundlogk[A_J]++;
					break;
				case 'T':
					backgroundlogk[T_J]++;
					break;
				case 'C':
					backgroundlogk[C_J]++;
					break;
				case 'G':
					backgroundlogk[G_J]++;
					break;
			}
		}

		setBackground(backgroundlogk[A_J],backgroundlogk[C_J],backgroundlogk[T_J],backgroundlogk[G_J]);
	}





	PSSM(int _length,int _psuedoCount=DEFEAULT_PSEUDO_COUNT,int _N_penalty=DEFAULT_N_PENALTY):M(NULL),length(_length),pseudoCount(_psuedoCount),N_penalty(_N_penalty)
	{
		M=new int*[length];
		for(int i=0;i<length;i++)
		{
			M[i]=new int[NT_ALPHABET_SIZE+1];
			for(int j=0;j<NT_ALPHABET_SIZE;j++)
			{
				M[i][j]=_psuedoCount;

			}

			M[i][NT_ALPHABET_SIZE]=_psuedoCount*NT_ALPHABET_SIZE;
		}

		numSeq=pseudoCount*NT_ALPHABET_SIZE;
	}





	void addmer(const string& nmer)
	{
		for(int i=0;i<length;i++)
		{
			char cur=nmer[i];
			switch(cur)
			{
			case 'A':
				M[i][A_J]++;
				M[i][NT_ALPHABET_SIZE]++;
				break;
			case 'C':
				M[i][C_J]++;
				M[i][NT_ALPHABET_SIZE]++;
				break;
			case 'T':
				M[i][T_J]++;
				M[i][NT_ALPHABET_SIZE]++;
				break;
			case 'G':
				M[i][G_J]++;
				M[i][NT_ALPHABET_SIZE]++;
				break;
			default:
				break;
			}
		}

		numSeq++;
	}

	void finalize()
	{
		//convert to log score
		for(int i=0;i<length;i++)
		{
			for(int j=0;j<NT_ALPHABET_SIZE;j++)
			{
				M[i][j]=LOGMUL*log(double(M[i][j])/M[i][NT_ALPHABET_SIZE]);
			}
		}
	}

	int score(const char* c,int p)
	{
		int thisScore=0;

		for(int i=0;i<length;i++)
		{
			char nt=c[p+i];
			switch(nt)
			{
			case 'A':
				thisScore+=M[i][A_J]-backgroundlogk[A_J];
				break;
			case 'T':
				thisScore+=M[i][T_J]-backgroundlogk[T_J];
				break;
			case 'C':
				thisScore+=M[i][C_J]-backgroundlogk[C_J];
				break;
			case 'G':
				thisScore+=M[i][G_J]-backgroundlogk[G_J];
				break;
			default:
				//anything else need penalty
				thisScore+=N_penalty;
			}
		}

		return thisScore;
	}

	class ScoreSeqPosMap: public multimap<int, int, greater<int> >
	{
	public:
		string seq;
		int lengthMotif;
		ScoreSeqPosMap(const string& _seq,int _lengthMotif):seq(_seq),lengthMotif(_lengthMotif){}

		void add(int score, int p)
		{
			//if(p>200)
			//	cerr<<"adding score "<<score<<" at "<<p<<endl;
			this->insert(value_type(score,p));
		}

		string getSequenceAt(int p)
		{
			//cerr<<p<<" "<<lengthMotif<<" "<<seq.length()<<endl;
			return seq.substr(p,lengthMotif);
		}

		void printScorePosSeq(ostream& os,int topN=10000000)
		{
			int printed=0;
			for(iterator i=this->begin();i!=this->end();i++)
			{
				os<<i->first<<"\t@"<<i->second<<"\t"<<this->getSequenceAt(i->second)<<endl;
				printed++;
				if(printed>=topN)
					break;
			}


		}

		void removeOverlappingMatches()
		{
			vector<iterator> toErase;

			bool* overlapMask=new bool[seq.length()];

			for(int i=0;i<seq.length();i++)
				overlapMask[i]=false;

			//sorted desending on score, remove lower scores
			for(iterator i=this->begin();i!=this->end();i++)
			{
				int pos=i->first;
				//check any nt of the region of this motif hit already covered
				bool previouslyOverlapped=false;
				for(int j=pos;j<pos+lengthMotif;j++)
				{
					if(overlapMask[j])
					{

						previouslyOverlapped=true;
						break;

					}
				}

				if (previouslyOverlapped)
				{
					//add this to destruction
					toErase.push_back(i);

				}
				else
				{
					//set overlap mask
					for(int j=pos;j<pos+lengthMotif;j++)
					{
						overlapMask[j]=true;
					}
				}
			}

			delete[] overlapMask;

			//now delete everything in toErase list
			for(vector<iterator>::iterator iI=toErase.begin();iI!=toErase.end();iI++)
			{
				this->erase(*iI);
			}
		}

	};

	void score(const string& seq,multiset<int,greater<int> > & scoreSet,ScoreSeqPosMap* pScoreSeqPosMap=NULL,bool dontAddNegativeScore=false)
	{
		int lseq=seq.length();
		const char*pset=seq.c_str();
		for(int p=0;p<lseq-length+1;p++)
		{
			int thisScore=score(pset,p);
			if(dontAddNegativeScore && thisScore<0)
			{
				continue;
			}
			scoreSet.insert(thisScore);
			if(pScoreSeqPosMap)
			{
				pScoreSeqPosMap->add(thisScore,p);
			}
		}
	}



	inline int MAX(int a,int b)
	{
		return (a>b)?a:b;
	}

	void print(ostream& os)
	{
		os<<"Model"<<endl;
		os<<"pos";
		for(int i=0;i<length;i++)
		{
			os<<"\t"<<i;
		}
		os<<endl;
		for(int j=0;j<NT_ALPHABET_SIZE;j++)
		{
			os<<ALPHABET_STRING[j];
			for(int i=0;i<length;i++)
			{
				os<<"\t"<<M[i][j];
			}
			os<<endl;
		}

		os<<"background:"<<endl;
		for(int j=0;j<NT_ALPHABET_SIZE;j++)
		{
			os<<ALPHABET_STRING[j]<<"\t"<<this->backgroundlogk[j]<<endl;
		}

	}

	void printScore(multiset<int,greater<int> > &scores,int topN=10000000)
	{
		int printed=0;
		for(multiset<int,greater<int> >::iterator i=scores.begin();i!=scores.end();i++)
		{
			cerr<<*i<<" ";
			printed++;
			if(printed>=topN)
				break;
		}

		cerr<<endl;
	}





	/*int getNumSigHits( string seq,double pvalueCutOff,int nShuffle=10)
	{

		multiset<int,greater<int> > fgScores;
		multiset<int,greater<int> > bgScores;



		this->setBackground(seq);

		this->print(cerr);

		multimap<int,string,greater<int> > fgScoreSeqMap;
		score(seq,fgScores,&fgScoreSeqMap);

		cerr<<"top 10 fgScores: "<<endl;
		this->printScoreSeq(fgScoreSeqMap,10);

		for(int s=0;s<nShuffle;s++)
		{

			//StringUtil::shuffleStringInPlace(seq);
			std::random_shuffle(seq.begin(),seq.end());
			score(seq,bgScores);
			cerr<<"shuffle "<<s<<" "<<seq.substr(0,10)<<" top 10 scores: ";
			this->printScore(bgScores,10);
		}

		//now nullScores and fgScores are null set and test set score sorted in descending order
		//int prevn=0;
		//int prepvalue=0;
		multiset<int,greater<int> >::iterator fgI=fgScores.begin();
		multiset<int,greater<int> >::iterator bgI=bgScores.begin();

		int nullSize=bgScores.size();
		int maxDiveInNull=nullSize*pvalueCutOff;  //automatically floored, which is good.
		cerr<<"nullSize="<<nullSize<<" pvalueCutOff="<<pvalueCutOff<<" maxDiveInNull="<<maxDiveInNull<<endl;

		this->printScore(bgScores,maxDiveInNull);

		int fgDive=0;
		int bgDive=0;

		while(true)
		{

			while(bgI!=bgScores.end() && *bgI>=*fgI)
			{
				bgDive++;
				bgI++;
				if (bgDive>=maxDiveInNull) //for zero pvalue as well
				{
					if (*bgI==*fgI)
						return MAX(0,fgDive-1);
					else
						return fgDive;
				}
			}



			while(fgI!=fgScores.end() && *fgI>*bgI)
			{
				fgDive++;
				fgI++;
			}

			//prevn=fgDive;

			if(fgI==fgScores.end())
				break;



		}

		return fgDive;

	}*/




	pair<int,ScoreSeqPosMap*> getNumSigHitsFDR( string seq,double fdrCutOff,bool noOverlap,bool keepFgScoreSeqPosMap=false,int nShuffle=10,bool dontAddNegativeScoreToFG=true)
	{

		multiset<int,greater<int> > fgScores;
		multiset<int,greater<int> > bgScores;

		//double final_pvalue=0.0;

		this->setBackground(seq);


		ScoreSeqPosMap*fgScoreSeqPosMap=new ScoreSeqPosMap(seq,this->length);




		score(seq,fgScores,fgScoreSeqPosMap,dontAddNegativeScoreToFG);
		if(noOverlap)
			fgScoreSeqPosMap->removeOverlappingMatches();

		//cerr<<"top 10 fgScores: "<<endl;
		//this->printScoreSeq(fgScoreSeqMap,10);

		for(int s=0;s<nShuffle;s++)
		{

			//StringUtil::shuffleStringInPlace(seq);
			std::random_shuffle(seq.begin(),seq.end());
			score(seq,bgScores);
			//cerr<<"shuffle "<<s<<" "<<seq.substr(0,10)<<" top 10 scores: ";
			//this->printScore(bgScores,10);
		}

		//now nullScores and fgScores are null set and test set score sorted in descending order
		//int prevn=0;
		//int prepvalue=0;
		ScoreSeqPosMap::iterator fgI=fgScoreSeqPosMap->begin();
		multiset<int,greater<int> >::iterator bgI=bgScores.begin();

		int nullSize=bgScores.size();
		//int maxDiveInNull=nullSize*pvalueCutOff;  //automatically floored, which is good.
		//cerr<<"nullSize="<<nullSize<<" pvalueCutOff="<<pvalueCutOff<<" maxDiveInNull="<<maxDiveInNull<<endl;



		int fgSize=fgScoreSeqPosMap->size();

		int fgDive=0;
		int bgDive=0;

		pair<double,int> prevValue(0.0,0);

		while(true)
		{

			while(bgI!=bgScores.end() && *bgI>=fgI->first)
			{
				bgDive++;
				bgI++;
				/*if (bgDive>=maxDiveInNull) //for zero pvalue as well
				{
					if (*bgI==*fgI)
						return MAX(0,fgDive-1);
					else
						return fgDive;
				}*/

			}



			while(fgI!=fgScoreSeqPosMap->end() && fgI->first>*bgI)
			{
				fgDive++;
				fgI++;
			}

			//prevn=fgDive;
			double pval=float(bgDive)/nullSize;
			double FDR=pval*fgSize/fgDive;

			//fdrfgDive.push_back(pair<double,int>(FDR,fgDive));

			if(FDR>=fdrCutOff)
			{
				fgDive= prevValue.second;
				break;
			}


			if(fgI==fgScoreSeqPosMap->end())
				break;

			prevValue.first=FDR;
			prevValue.second=fgDive;

		}

		//cerr<<"final FDR="<<prevValue.first;
		//cerr<<" final pvalue="<<(prevValue.first*fgDive/fgSize)<<endl;
		//this->printScore(bgScores,bgDive);

		if(keepFgScoreSeqPosMap)
			return pair<int,ScoreSeqPosMap*>(fgDive,fgScoreSeqPosMap);


		delete fgScoreSeqPosMap;

		return pair<int,ScoreSeqPosMap*>(fgDive,NULL);

	}

	/*int getNumSigHitsBenCorr(const string &seq, double pvalueCutOff,int nShuffle=10)
	{
		return getNumSigHits(seq,pvalueCutOff/seq.length(),nShuffle);
	}*/

	~PSSM(){

		if(M)
		{

			for(int i=0;i<length;i++)
				delete[] M[i];

			delete[] M;
		}

	}

};

class MotifFamily
{

public:
	map<string,vector<string>* > NmerList;
	map<string,vector<string>* > NmerAlign;
	int maxLenAlign;
	int maxLenNoGap;

	//remember to free the PSSM after use!
	PSSM* getPSSMOf(const string& familyName)
	{
		map<string,vector<string>* >::iterator i=NmerAlign.find(familyName);
		if (i==NmerAlign.end())
			return NULL;

		vector<string>* nmers=i->second;

		if(nmers->size()==0)
			return NULL;

		PSSM* pssm=new PSSM(nmers->at(0).length());

		for(vector<string>::iterator j=nmers->begin();j!=nmers->end();j++)
		{
			pssm->addmer(*j);
		}
		pssm->finalize();

		return pssm;
	}

	class SortDescendingLength

	{
	public:
		bool operator()(const string& i,const string& j)
		{
			return i.length()>j.length();
		}

	} sortDescendingLength;

	void finalize()
	{
		for(map<string,vector<string>* >::iterator i=NmerList.begin();i!=NmerList.end();i++)
		{
			vector<string>* vstr=i->second;
			sort(vstr->begin(),vstr->end(),sortDescendingLength);
		}
	}



	MotifFamily(string filename,int NmerOnly=-1)
	{
		//string curName=""
		ifstream fin(filename.c_str());
		char *line=new char[_BUFFER_SIZE];

		vector<string> splits;
		line[0]='\0';
		maxLenAlign=0;
		maxLenNoGap=0;

		vector<string> splits2;

		while(true)
		{
			fin.getline(line,_BUFFER_SIZE);

			//cerr<<"line:"<<line<<endl;

			if(line[0]=='\0')
				break;

			StringUtil::splitNoEmpty(line," ",splits);
			if(splits.size()<2)
				break;

			string sequence=splits[1];
			StringUtil::splitNoEmpty(splits[0],"_",splits2);
			string familyName=splits2[0];

			string nogapsequence=StringUtil::replace(splits[1],"-","");
			if(NmerOnly!=-1 && NmerOnly!=nogapsequence.length())
			{
				continue;
			}

			vector<string>* curList;
			vector<string>* curAlign;

			map<string,vector<string>* >::iterator f=NmerList.find(familyName);
			if(f==NmerList.end())
			{
				curList=new vector<string>;
				NmerList.insert(map<string,vector<string>* >::value_type(familyName,curList) );
			}else
				curList=f->second;

			f=NmerAlign.find(familyName);
			if(f==NmerAlign.end())
			{
				curAlign=new vector<string>;
				NmerAlign.insert(map<string,vector<string>* >::value_type(familyName,curAlign));
			}else
				curAlign=f->second;

			if(nogapsequence.length()>maxLenNoGap)
				maxLenNoGap=nogapsequence.length();
			if(sequence.length()>maxLenAlign)
				maxLenAlign=sequence.length();

			curList->push_back(nogapsequence);
			curAlign->push_back(sequence);



		}

		delete[] line;
		finalize();
	}



	~MotifFamily()
	{
		for(map<string,vector<string>* >::iterator i=NmerList.begin();i!=NmerList.end();i++)
		{
			delete i->second;
		}
		for(map<string,vector<string>* >::iterator i=NmerAlign.begin();i!=NmerAlign.end();i++)
		{
			delete i->second;
		}
	}

};

class PSSMSet: public map<string,PSSM*> {
public:

	PSSMSet(MotifFamily& motifFamily)
	{
		for(map<string,vector<string>* >::iterator i=motifFamily.NmerAlign.begin();i!=motifFamily.NmerAlign.end();i++)
		{
			string famName=i->first;
			//vector<string>* aligns=i->second;
			PSSM* pssm=motifFamily.getPSSMOf(famName);
			this->insert(map<string,PSSM*>::value_type(famName,pssm));

		}
	}

	~PSSMSet()
	{
		for(map<string,PSSM*>::iterator i=this->begin();i!=this->end();i++)
		{
			delete i->second;
		}
	}

	/*int getNumSigHitsFDR(const string& seq,map<string,PSSM*>::iterator pssmI,double FDRCutOff,int nShuffle=100)
	{
		PSSM* pssm=pssmI->second;
		pssm->setBackground(seq);
		return pssm->getNumSigHitsFDR(seq,FDRCutOff,nShuffle);

	}
	int getNumSigHitsFDR(const string&seq, const string& familyName,double FDRCutOff,int nShuffle=100)
	{
		map<string,PSSM*>::iterator pssmI=this->find(familyName);
		if(pssmI==this->end())
			return -1;

		return this->getNumSigHitsFDR(seq,pssmI,FDRCutOff,nShuffle);
	}*/
};


void JHypergeometric_Pvalue(const vector<int> &input,vector<double>& result,double ErrValue=-1)
{
	char tmpfilename[100];
	char floatResult[100];
	tmpnam(tmpfilename);
	cerr<<"tmpfilename="<<tmpfilename<<endl;
	ofstream fout(tmpfilename);
	for(vector<int>::const_iterator i=input.begin();i!=input.end();)
	{
		//samt sam popt popt
		int samt=*(i++);
		int sam=*(i++);
		int popt=*(i++);
		int pop=*(i++);
		fout<<StringUtil::str(samt)<<" "<<StringUtil::str(sam)<<" "<<StringUtil::str(popt)<<" "<<StringUtil::str(pop)<<endl;
		if(i==input.end())
			break;
	}
	fout.close();

	//now call Java Helper

	string command=string("java HypergeometricPvalue ")+tmpfilename;

	FILE *fp;
	int status;



	fp = popen(command.c_str(), "r");
	if (fp == NULL){
		cerr<<"Cannot open Java helper program "<<command<<endl;
		exit(0);
	}

	double fresult;

	while (fgets(floatResult, 100, fp) != NULL){
		if(floatResult[0]=='_')
			fresult=ErrValue;
		else
			fresult=StringUtil::atof(floatResult);

		result.push_back(fresult);
	}

	status = pclose(fp);


	remove(tmpfilename);


}

void putInputVector(vector<int>& input, int samt,int sam, int popt, int pop)
{
	input.push_back(samt);
	input.push_back(sam);
	input.push_back(popt);
	input.push_back(pop);
}

#endif

class WordFreqMap: public map<string, int>
{
	
};

class Bin
{

public :
 vector<string> sequences;
};

class WordResult
{
public:
	
	double pvalue;
	double expectedProb;
	int expected;
	double observedProb;
	int observed;
	string word;
	
	inline WordResult(double _pvalue,double _expectedProb,int _expected,double _observedProb,int _observed,string _word):
		pvalue(_pvalue),expectedProb(_expectedProb),expected(_expected),observedProb(_observedProb),observed(_observed),word(_word){}
};

inline string convertT2U(const string& from)
{
	string result(from);
	int l=result.length();
	for(int i=0;i<l;i++)
	{
		if(result[i]=='T')
			result[i]='U';
	}
	
	return result;
}


class WordArray:public vector<WordResult>
{
	
};

class PvalueWordMap: public map<double,WordArray>
{
	
};

class Bins: public map<int,Bin>
{
	
};

class ElementStruct
{
public:
	WordFreqMap wordreq;
	int countNmers;
	Bins bins;
	
	inline ElementStruct(): countNmers(0){}
};



class HeaderStruct: public map<string, ElementStruct>
{
	
};

class DictNmer: public map<string, HeaderStruct>
{
	
};

double CGValue(const string&seq)
{
	int CG=0;
	int lreal=0;
	int lseq=seq.length();
	for(int i=0;i<lseq;i++)
	{
		char c=seq[i];
		if(c=='C' || c=='G')
				CG++;
		else if(c!='A' and c!='T')
			 continue;
		
		lreal++;
	}
	
	return float(CG)/lreal;
}




inline double min(double x, double y)
{
	return (x<y)?x:y;
}


inline int randrange(int from0,int to1)
{
	return rand()%(to1-from0)+from0;
}

bool hasN(const string& s)
{
	const char* c_str=s.c_str();
	while(*c_str!='\0')
	{
		char c=*(c_str++);
		if(c!='A' && c!='G' && c!='T' && c!='C')
			return true;
		
	}
	
	return false;
}

void hashSequenceSimple(ElementStruct& element,const string& content,int N,bool noOverlap)
{
	int	lcontent=content.length();
	int end1=lcontent-N+1;
	
	int countN=0;
	int countNmers=0;
	
	WordFreqMap& wordFreq=element.wordreq;
	
	map<string,int> block;


	for(int i=0;i<end1;i++)
	{
		string wordToHash=content.substr(i,N);
		if(hasN(wordToHash))
			continue;
		
		map<string,int>::iterator blockI;

		if(noOverlap)
		{
			blockI=block.find(wordToHash);
			if(blockI==block.end())
			{
				pair<map<string,int>::iterator,bool > result=block.insert(map<string,int>::value_type(wordToHash,0));
				blockI=result.first;
			}

			if(i<blockI->second)
			{
				continue; //blocked by overlap
			}
		}

		if(!noOverlap)
			countNmers++;
		
		WordFreqMap::iterator WordFreqI=wordFreq.find(wordToHash);
		if(WordFreqI==wordFreq.end())
		{
			pair<WordFreqMap::iterator,bool> result=wordFreq.insert(WordFreqMap::value_type(wordToHash,0));
			WordFreqI=result.first;

		}
		
		
		WordFreqI->second++;

		if(noOverlap)
			blockI->second=i+N;


	}

	if(noOverlap)
	{
		countNmers=lcontent/N;
	}
	
	//cerr<<"B"<<endl;
	element.countNmers+=countNmers;
	
}

void hashSequenceSimpleForPSSM(ElementStruct& element,const string& content,int N,MotifFamily *motifFam,bool noOverlap,int PSSMnShuffle,double PSSMFDRThreshold)
{
	int	lcontent=content.length();
	int end1=lcontent-N+1;

	int countN=0;
	int countNmers=0;

	WordFreqMap& wordFreq=element.wordreq;

	//maxN=0;


	int NmerListSize=motifFam->NmerList.size();

	int* blocker=new int[NmerListSize];

	for(int i=0;i<NmerListSize;i++)
		blocker[i]=0;


	if(!noOverlap){
		for(int i=0;i<end1;i++)
		{
			string wordToHash=content.substr(i,N);
			if(hasN(wordToHash))
				continue;


			countNmers++;
		}
	}

	for(map<string,vector<string>* >::iterator j=motifFam->NmerList.begin();j!=motifFam->NmerList.end();j++)
	{


		string familyName=j->first;

		PSSM* pssm=motifFam->getPSSMOf(familyName);

		if(!pssm)
		{
			cerr<<"PSSM creation error"<<endl;
			continue;
		}

		pssm->setBackground(content);



		pair<int,PSSM::ScoreSeqPosMap*> searchResult=pssm->getNumSigHitsFDR(content,PSSMFDRThreshold,true,true,PSSMnShuffle);
		//int numSig=pssm->getNumSigHits(seq,0.05,100);
		int numSig=searchResult.first;
		PSSM::ScoreSeqPosMap *hits=searchResult.second;

		//cerr<<"numSig="<<numSig<<endl;
		//hits->printScorePosSeq(cerr,numSig);

		delete hits;



		if(pssm)
			delete pssm;

		WordFreqMap::iterator WordFreqI=wordFreq.find(familyName);
		if(WordFreqI==wordFreq.end())
		{
			wordFreq.insert(WordFreqMap::value_type(familyName,numSig));

		}else
		{
			//cerr<<"repeat of family??"<<endl;
			//exit(1);
			WordFreqI->second+=numSig;
		}



	}








	if(noOverlap)
		countNmers=lcontent/motifFam->maxLenAlign;
	//cerr<<"B"<<endl;
	element.countNmers+=countNmers;

}

void hashSequenceSimpleForFam(ElementStruct& element,const string& content,int N,MotifFamily *motifFam,bool noOverlap)
{
	int	lcontent=content.length();
	int end1=lcontent-N+1;

	int countN=0;
	int countNmers=0;

	WordFreqMap& wordFreq=element.wordreq;


	int NmerListSize=motifFam->NmerList.size();

	int* blocker=new int[NmerListSize];

	for(int i=0;i<NmerListSize;i++)
		blocker[i]=0;

	for(int i=0;i<end1;i++)
	{
		string wordToHash=content.substr(i,N);
		if(hasN(wordToHash))
			continue;

		if(!noOverlap)
			countNmers++;

		int blockerI=-1;
		for(map<string,vector<string>* >::iterator j=motifFam->NmerList.begin();j!=motifFam->NmerList.end();j++)
		{

			blockerI++;

			if (noOverlap && i<blocker[blockerI])
			{
				//for this family of nmers, already mapped and blocked for overlap
				continue;
			}

			string familyName=j->first;
			vector<string>* nmers=j->second;

			WordFreqMap::iterator WordFreqI=wordFreq.find(familyName);
			if(WordFreqI==wordFreq.end())
			{
				pair<WordFreqMap::iterator, bool > result=wordFreq.insert(WordFreqMap::value_type(familyName,0));
				WordFreqI=result.first;
			}

			for(vector<string>::iterator k=nmers->begin();k!=nmers->end();k++)
			{
				if (*k==wordToHash)
				{
					WordFreqI->second++;
					if(noOverlap){
						blocker[blockerI]=i+k->length();
						break;
					}
				}
			}



		}



	}




	if(noOverlap)
		countNmers=lcontent/N;
	//cerr<<"B"<<endl;
	element.countNmers+=countNmers;

}

/*void hashSequenceSimpleForFamPSSM(ElementStruct& element,const string& content,int N,PSSMSet *pssmset,bool noOverlap)
{
	int	lcontent=content.length();
	int end1=lcontent-N+1;

	int countN=0;
	int countNmers=0;

	WordFreqMap& wordFreq=element.wordreq;


	int NmerListSize=motifFam->NmerList.size();

	int* blocker=new int[NmerListSize];

	for(int i=0;i<NmerListSize;i++)
		blocker[i]=0;

	for(int i=0;i<end1;i++)
	{
		string wordToHash=content.substr(i,N);
		if(hasN(wordToHash))
			continue;

		if(!noOverlap)
			countNmers++;

		int blockerI=-1;
		for(map<string,vector<string>* >::iterator j=motifFam->NmerList.begin();j!=motifFam->NmerList.end();j++)
		{

			blockerI++;

			if (noOverlap && i<blocker[blockerI])
			{
				//for this family of nmers, already mapped and blocked for overlap
				continue;
			}

			string familyName=j->first;
			vector<string>* nmers=j->second;

			WordFreqMap::iterator WordFreqI=wordFreq.find(familyName);
			if(WordFreqI==wordFreq.end())
			{
				pair<WordFreqMap::iterator, bool > result=wordFreq.insert(WordFreqMap::value_type(familyName,0));
				WordFreqI=result.first;
			}

			for(vector<string>::iterator k=nmers->begin();k!=nmers->end();k++)
			{
				if (*k==wordToHash)
				{
					WordFreqI->second++;
					if(noOverlap){
						blocker[blockerI]=i+k->length();
						break;
					}
				}
			}



		}



	}




	if(noOverlap)
		countNmers=lcontent/N;
	//cerr<<"B"<<endl;
	element.countNmers+=countNmers;

}*/




void constructHGWEProportionBGPerElement(ElementStruct& FGElement,ElementStruct& BGpElement,ElementStruct& BGElement,int N,MotifFamily*motifFamily,bool noOverlap,bool usePSSM,int PSSMnShuffle,double PSSMFDRThreshold)
{
	Bins& FGBins=FGElement.bins;
	Bins& BGpBins=BGpElement.bins;
	
	double minRatio=10000000;
	
	for(Bins::iterator FGBinI=FGBins.begin();FGBinI!=FGBins.end();FGBinI++)
	{
		int CGBinKey=FGBinI->first;
		Bin& CGBinFG=FGBinI->second;
		Bins::iterator BGpBinI=BGpBins.find(CGBinKey);
		
		cerr<<"getting bg/fg ratio for bin "<<CGBinKey<<" ";
		
		if(BGpBinI==BGpBins.end())
		{
			cerr<<"that bin is not found in bacakground"<<endl;
		}
		else
		{
			Bin& CGBinBGp=BGpBinI->second;
			
			vector<string>& FgSequences=CGBinFG.sequences;
			int nCGBinFg=FgSequences.size();
			int nCGBinBgp=CGBinBGp.sequences.size();
			
			for(vector<string>::iterator i=FgSequences.begin();i!=FgSequences.end();i++)
			{
				//cerr<<"hashing sequence "<<*i<<endl;

				if(motifFamily){
					if(usePSSM)
					{
						hashSequenceSimpleForPSSM(FGElement,*i,N,motifFamily,noOverlap,PSSMnShuffle,PSSMFDRThreshold);
						hashSequenceSimpleForPSSM(BGElement,*i,N,motifFamily,noOverlap,PSSMnShuffle,PSSMFDRThreshold); //the same as direct copy!!
					}else
					{
						hashSequenceSimpleForFam(FGElement,*i,N,motifFamily,noOverlap);
						hashSequenceSimpleForFam(BGElement,*i,N,motifFamily,noOverlap); //the same as direct copy!!
					}
				}else
				{
					hashSequenceSimple(FGElement,*i,N,noOverlap);
					hashSequenceSimple(BGElement,*i,N,noOverlap);//the same as direct copy!!
				}
			}
			double thisRatio=double(nCGBinBgp)/double(nCGBinFg);
			minRatio=min(minRatio,thisRatio);
			cerr<< thisRatio<<endl;
		}
		
	}
	
	cerr<<"minRatio="<<minRatio<<endl;
	
	for(Bins::iterator FGBinI=FGBins.begin();FGBinI!=FGBins.end();FGBinI++)
	{
		int CGBinKey=FGBinI->first;
		Bin& CGBinFG=FGBinI->second;
		Bins::iterator BGpBinI=BGpBins.find(CGBinKey);
		
		//cerr<<"getting bg/fg ratio for bin "<<CGBinKey<<" ";
		
		if(BGpBinI!=BGpBins.end())
		{
			Bin& CGBinBGp=BGpBinI->second;
			vector<string>& BgpSequences=CGBinBGp.sequences;
			int nCGBinBgp=BgpSequences.size();
			int nCGBinFG=CGBinFG.sequences.size();
			
			int requiredNo=int(minRatio*nCGBinFG);
			if(requiredNo>nCGBinBgp)
			{
				cerr<<"strange Error: requiredNo > nCGBinBgp"<<endl;
			}
			else if(requiredNo==nCGBinBgp)
			{
				for(vector<string>::iterator i=BgpSequences.begin();i!=BgpSequences.end();i++)
				{
					if(motifFamily)
					{
						if(usePSSM)
						{
							hashSequenceSimpleForPSSM(BGElement,*i,N,motifFamily,noOverlap,PSSMnShuffle,PSSMFDRThreshold);
						}else
						{


							hashSequenceSimpleForFam(BGElement,*i,N,motifFamily,noOverlap);
						}
					}else
					{
						hashSequenceSimple(BGElement,*i,N,noOverlap);
					}

				}
			}
			else
			{
				for(int i=0;i<requiredNo;i++)
				{
					int randi=randrange(0,nCGBinBgp-i);

					if(motifFamily)
					{

						if(usePSSM){
							hashSequenceSimpleForPSSM(BGElement,BgpSequences[randi],N,motifFamily,noOverlap,PSSMnShuffle,PSSMFDRThreshold);

						}else
						{
							hashSequenceSimpleForFam(BGElement,BgpSequences[randi],N,motifFamily,noOverlap);
						}
					}
					else
					{
						hashSequenceSimple(BGElement,BgpSequences[randi],N,noOverlap);
					}
					BgpSequences[randi]=BgpSequences[nCGBinBgp-i-1];
				}
			}
		}
	}
	
}

void constructHGWEProportionBg(DictNmer& FG,DictNmer& BGp,DictNmer& BG,int N,MotifFamily *motifFamily,bool noOverlap,bool usePSSM, int PSSMnShuffle, double PSSMFDRThreshold)
{
	DictNmer::iterator HeaderStructI;
	
	for(HeaderStructI=FG.begin();HeaderStructI!=FG.end();HeaderStructI++)
	{
		string header=HeaderStructI->first;
		HeaderStruct& perHeaderFG=HeaderStructI->second;
		BG.insert(DictNmer::value_type(header,HeaderStruct()));
		
		HeaderStruct& perHeaderBGp=BGp.find(header)->second;
		HeaderStruct& perHeaderBG=BG.find(header)->second;
		
		HeaderStruct::iterator ElementStructI;
		for(ElementStructI=perHeaderFG.begin();ElementStructI!=perHeaderFG.end();ElementStructI++)
		{
			string elementKey=ElementStructI->first;
			ElementStruct& perElementFG=ElementStructI->second;
			
			HeaderStruct::iterator ElementStructIBGp=perHeaderBGp.find(elementKey);
			ElementStruct& perElementBGp=ElementStructIBGp->second;
			
			
			HeaderStruct::iterator ElementStructIBG=perHeaderBG.find(elementKey);
			if(ElementStructIBG==perHeaderBG.end())
			{
				pair<HeaderStruct::iterator, bool> insertStat=perHeaderBG.insert(HeaderStruct::value_type(elementKey,ElementStruct()));
				ElementStructIBG=insertStat.first;
			}
			ElementStruct& perElementBG=ElementStructIBG->second;
			
			
			cerr<<"construct Proportionate BG for element:"<<header<<":"<<elementKey<<":"<<endl;
			constructHGWEProportionBGPerElement(perElementFG,perElementBGp,perElementBG,N,motifFamily,noOverlap, usePSSM, PSSMnShuffle, PSSMFDRThreshold);
			
			
		}
	}
}


void hashSequenceToCGBinByIntervals(DictNmer& dict,string filename,int headerRow1,int startRow1,int strandCol1,int fromCol1,int toCol1,double CGBinInterval)
{
	fstream fin(filename.c_str());
	
	int fromCol0=fromCol1-1;
	int strandCol0=strandCol1-1;
	
	vector<string> header;
	
	int lino=0;
	
	cerr<<"hashing stuffaaa"<<endl;
	
	char *line=new char[_BUFFER_SIZE];
	
	line[0]='\0';
	
	
	
	vector<string> spliton;
	
	while(true)
	{
		fin.getline(line,_BUFFER_SIZE);
		
		//cerr<<"line:"<<line<<endl;
		
		if(line[0]=='\0')
			break;
		
		lino++;
		StringUtil::split(line,"\t",spliton);
		if(lino==headerRow1)
		{
			for(int i=fromCol0;i<toCol1;i++)
			{
				header.push_back(spliton[i]);
				cerr<<"New Header "<<spliton[i]<<endl;
				dict.insert(DictNmer::value_type(spliton[i],HeaderStruct()));
			}
		}
		
		if(lino<startRow1)
		
		
			continue;
		
	
	
	for(int i=fromCol0;i<toCol1;i++)
	{
		string key=header[i-fromCol0];
		string strand=spliton[strandCol0];
		string content=StringUtil::toUpper(spliton[i]);
		if(content.length()<3)
			continue;
		
		vector<string> contentSplits;
		StringUtil::split(content,"|",contentSplits);
		
		DictNmer::iterator HeaderStructI=dict.find(key);
		HeaderStruct& perHeader=HeaderStructI->second;
		
		
		for(vector<string>::iterator j=contentSplits.begin();j!=contentSplits.end();j++)
		{
			string contentSpliton=*j;
			if(contentSpliton.length()<3)
				continue;
			
			string elementKey=contentSpliton.substr(0,2);
			string elementSeq=contentSpliton.substr(2);
			
			
			
			//if(strand=="-")
			//	elementSeq=reverse_complement(elementSeq);
			
			HeaderStruct::iterator ElementStructI=perHeader.find(elementKey);
			
			if(ElementStructI==perHeader.end())
			{
				cerr<<"new element:"<<key<<":"<<elementKey<<endl;
				perHeader.insert(HeaderStruct::value_type(elementKey,ElementStruct()));
				ElementStructI=perHeader.find(elementKey);
			}
			
			ElementStruct& perElement=ElementStructI->second;
			
			Bins& CGBins=perElement.bins;
			double CGVal=CGValue(elementSeq);
			int CGBinKey=int(CGVal/CGBinInterval);
			
			Bins::iterator BinI=CGBins.find(CGBinKey);
			
			if(BinI==CGBins.end())
			{
				CGBins.insert(Bins::value_type(CGBinKey,Bin()));
				BinI=CGBins.find(CGBinKey);
			}
			
			Bin& bin=BinI->second;
			
			bin.sequences.push_back(elementSeq);
			
			
			
		}
	}
	
	
	
	}
	
	
	delete[] line;
	fin.close();
}




void findEnrichedMotifByHGWE(string filenamefg,string filenamebg,int headerRow1,int startRow1,int strandCol1,int fromCol1,int toCol1,int N,double CGBinInterval,double FDRCutOff,double pvalueCutOff,int topCount,string FamilyFile,bool noOverlap,bool usePSSM,int PSSMnShuffle,double PSSMFDRThreshold)
{

	MotifFamily *motifFam=NULL;



	if(FamilyFile!=".")
	{
		cerr<<">>motif family mode"<<endl;
		motifFam=new MotifFamily(FamilyFile,N);
		if(usePSSM)
		{
			cerr<<">>PSSM mode"<<endl;
		}
	}
	else
	{
		if(usePSSM) //but no family file
		{
			cerr<<"use PSSM but no family file defined. Abort";
			return;
		}
	}

	if(noOverlap)
	{
		cerr<<"no overlap mode"<<endl;
	}


	DictNmer DictNmerFG;
	DictNmer DictNmerBGp;
	DictNmer DictNmerBG;
	
	cerr<<"put sequence into CG value bins"<<endl;
	
	cout <<"Header\tElement\tp-value\tFDR\texpected rate\texpected freq\tobserved rate\tobserved freq\tword"<<endl;

	
	hashSequenceToCGBinByIntervals(DictNmerFG,filenamefg,headerRow1,startRow1,strandCol1,fromCol1,toCol1,CGBinInterval);
	hashSequenceToCGBinByIntervals(DictNmerBGp,filenamebg,headerRow1,startRow1,strandCol1,fromCol1,toCol1,CGBinInterval);	
	
	
	cerr<<"construct appropriate BG"<<endl;
	
	constructHGWEProportionBg(DictNmerFG,DictNmerBGp,DictNmerBG,N,motifFam,noOverlap,usePSSM,PSSMnShuffle,PSSMFDRThreshold);
	
	//now the real stuff
	
	DictNmer::iterator HeaderStructI;
	
	for(HeaderStructI=DictNmerFG.begin();HeaderStructI!=DictNmerFG.end();HeaderStructI++)
	{
		string header=HeaderStructI->first;
		
		HeaderStruct& perHeaderFG=HeaderStructI->second;
		//HeaderStruct& perHeaderBGp=BGp.find(header)->second;
		HeaderStruct& perHeaderBG=DictNmerBG.find(header)->second;
		
		cerr<< header<<":"<<endl;
		
		HeaderStruct::iterator ElementStructI;
		


		for(ElementStructI=perHeaderFG.begin();ElementStructI!=perHeaderFG.end();ElementStructI++)
		{
			cerr<<"a"<<endl;
			string elementKey=ElementStructI->first;
			ElementStruct& perElementFG=ElementStructI->second;
			ElementStruct& perElementBG=perHeaderBG.find(elementKey)->second;
			cerr<< header<<":"<<elementKey<<":"<<endl;
			//cout<<header<<":"<<elementKey<<":"<<endl;
			WordFreqMap& wordFreqMapFG=perElementFG.wordreq;
			WordFreqMap& wordFreqMapBG=perElementBG.wordreq;
			cerr<<"b"<<endl;
			int numWordsFG=wordFreqMapFG.size();
			int numWordsBG=wordFreqMapBG.size();
			int numNmersInFG=perElementFG.countNmers;
			int numNmersInBG=perElementBG.countNmers;
			PvalueWordMap pvalueWordMap;
			cerr<<"c"<<endl;
			
#ifdef USE_JAVA_HYPERGEOMETRIC
			vector<double> pvalues;
			vector<int> jInputVector;


			/*if(motifFam)
			{


				for(map<string,vector<string>* >::iterator famI=motifFam->NmerList.begin();famI!=motifFam->NmerList.end();famI++)
				{
					string familyName=famI->first;
					vector<string>* nmers=famI->second;
					int freqFG=0;
					int freqBG=0;

					for(vector<string>::iterator nmerI=nmers->begin();nmerI!=nmers->end();nmerI++)
					{
						string word=*nmerI;
						WordFreqMap::iterator wordFreqIFG=wordFreqMapFG.find(word);
						WordFreqMap::iterator wordFreqIBG=wordFreqMapBG.find(word);
						if(wordFreqIFG==wordFreqMapFG.end() || wordFreqIBG==wordFreqMapBG.end())
						{
							continue;
						}

						freqFG+=wordFreqIFG->second;
						freqBG+=wordFreqIBG->second;





					::putInputVector(jInputVector,freqFG,numNmersInFG,freqBG,numNmersInBG);
				}
			}
			else
			{*/
				for(WordFreqMap::iterator WordIFG=wordFreqMapFG.begin();WordIFG!=wordFreqMapFG.end();WordIFG++)
				{
					cerr<<"d"<<endl;
					string word=WordIFG->first;
					int freqFG=WordIFG->second;
					int freqBG=wordFreqMapBG.find(word)->second;
					cerr<<"e"<<endl;
					//cerr<< word <<" : Hypergeometric_cdf(m,mt,n,nt): "<<" "<<numNmersInBG<<" " <<freqBG<< " "<<numNmersInFG<<" " << freqFG<< " ";
					//double pvalue=pvalue_enrichment(numNmersInBG,freqBG,numNmersInFG,freqFG);
					//cerr<<" pvalue="<<pvalue<<endl;
					::putInputVector(jInputVector,freqFG,numNmersInFG,freqBG,numNmersInBG);
					cerr<<"f"<<endl;
				}
			//}

			//now calculate
				cerr<<"g"<<endl;
			::JHypergeometric_Pvalue(jInputVector,pvalues,1);
			cerr<<"h"<<endl;
			vector<double>::iterator pvaluesI=pvalues.begin();

			//second pass

			/*if(motifFam)
			{
				int jInputVectorI=0;

				for(map<string,vector<string>* >::iterator famI=motifFam->NmerList.begin();famI!=motifFam->NmerList.end();famI++)
				{
					string familyName=famI->first;
					vector<string>* nmers=famI->second;

					int freqFG=jInputVector[jInputVectorI];

					int freqBG=jInputVector[jInputVectorI+2];


					//::putInputVector(jInputVector,freqFG,numNmersInFG,freqBG,numNmersInBG);

					jInputVectorI+=4;

					cerr<< familyName <<" : Hypergeometric_cdf(m,mt,n,nt): "<<" "<<numNmersInBG<<" " <<freqBG<< " "<<numNmersInFG<<" " << freqFG<< " ";
					double pvalue=(*pvaluesI++);
					cerr<<" pvalue="<<pvalue<<endl;

					PvalueWordMap::iterator wordArrayI=pvalueWordMap.find(pvalue);

					if(wordArrayI==pvalueWordMap.end())
					{
						pair<PvalueWordMap::iterator,bool> insertStat=pvalueWordMap.insert(PvalueWordMap::value_type(pvalue,WordArray()));
						wordArrayI=insertStat.first;
					}

					WordArray& wordArray=wordArrayI->second;

					wordArray.push_back(WordResult(pvalue,double(freqBG)/numNmersInBG,freqBG,double(freqFG)/numNmersInFG,freqFG,familyName));

				}
			}
			else
			{*/
				for(WordFreqMap::iterator WordIFG=wordFreqMapFG.begin();WordIFG!=wordFreqMapFG.end();WordIFG++)
				{
					cerr<<"i"<<endl;
					string word=WordIFG->first;
					int freqFG=WordIFG->second;
					int freqBG=wordFreqMapBG.find(word)->second;

					cerr<< word <<" : Hypergeometric_cdf(m,mt,n,nt): "<<" "<<numNmersInBG<<" " <<freqBG<< " "<<numNmersInFG<<" " << freqFG<< " ";
					double pvalue=(*pvaluesI++);
					cerr<<" pvalue="<<pvalue<<endl;

					PvalueWordMap::iterator wordArrayI=pvalueWordMap.find(pvalue);

					if(wordArrayI==pvalueWordMap.end())
					{
						pair<PvalueWordMap::iterator,bool> insertStat=pvalueWordMap.insert(PvalueWordMap::value_type(pvalue,WordArray()));
						wordArrayI=insertStat.first;
					}

					cerr<<"j"<<endl;
					WordArray& wordArray=wordArrayI->second;
					cerr<<"k"<<endl;
					wordArray.push_back(WordResult(pvalue,double(freqBG)/numNmersInBG,/*freqBG*/double(freqBG)/numNmersInBG*numNmersInFG,double(freqFG)/numNmersInFG,freqFG,word));
					cerr<<"l"<<endl;
				}
			//}


#else
			
			
			for(WordFreqMap::iterator WordIFG=wordFreqMapFG.begin();WordIFG!=wordFreqMapFG.end();WordIFG++)
			{
				string word=WordIFG->first;
				int freqFG=WordIFG->second;
				int freqBG=wordFreqMapBG.find(word)->second;
				
				cerr<< word <<" : Hypergeometric_cdf(m,mt,n,nt): "<<" "<<numNmersInBG<<" " <<freqBG<< " "<<numNmersInFG<<" " << freqFG<< " ";
				double pvalue=pvalue_enrichment(numNmersInBG,freqBG,numNmersInFG,freqFG);
				cerr<<" pvalue="<<pvalue<<endl;
				
				PvalueWordMap::iterator wordArrayI=pvalueWordMap.find(pvalue);
				
				if(wordArrayI==pvalueWordMap.end())
				{
					pair<PvalueWordMap::iterator,bool> insertStat=pvalueWordMap.insert(PvalueWordMap::value_type(pvalue,WordArray()));
					wordArrayI=insertStat.first;
				}
				
				WordArray& wordArray=wordArrayI->second;
				
				wordArray.push_back(WordResult(pvalue,double(freqBG)/numNmersInBG,freqBG,double(freqFG)/numNmersInFG,freqFG,word));
			}	
			
#endif

			
				//cout<<"#NmersFG: "<< numNmersInFG<<" #NmersBG:" <<numNmersInBG<<" #wordsFG:"<< numWordsFG <<" #wordsBG:"<< numWordsBG<<endl;
				int outputed=0;
				
				int familySize=0;
				if(motifFam)
					familySize=motifFam->NmerList.size();


				for(PvalueWordMap::iterator wordArrayI=pvalueWordMap.begin();wordArrayI!=pvalueWordMap.end();wordArrayI++)
				{
					cerr<<"m"<<endl;
					double pvalue=wordArrayI->first;
					WordArray& wordArray=wordArrayI->second;
					cerr<<"n"<<endl;
					if(pvalueCutOff>=0 && pvalue>pvalueCutOff)
						break;
					cerr<<"o"<<endl;
					int nWordWithThisPvalue=wordArray.size();
					if(outputed+nWordWithThisPvalue>topCount && outputed>0)
					{
						cerr<<"top Count cut-offed_:"<<(outputed+nWordWithThisPvalue)<<endl;
						break;
					}
					double FDR=0;
					
					for(WordArray::iterator WordResultI=wordArray.begin();WordResultI!=wordArray.end();WordResultI++)
					{
						cerr<<"p"<<endl;
						WordResult& wordResult=*WordResultI;
						if(motifFam)
							FDR=wordResult.pvalue*familySize/(outputed+1);
						else
							FDR=wordResult.pvalue*numWordsFG/(outputed+1);
						cerr<<"q"<<endl;

						if(FDR>1)
						 {
							 FDR=1.0;
						 }
						if(FDRCutOff>=0 && FDR>FDRCutOff)
							break;
						
						string outputWord;
						if(motifFam)
							outputWord=wordResult.word;
						else
							outputWord=convertT2U(wordResult.word);
						cerr<<"s"<<endl;
						cout<<header<<"\t"<<elementKey<<"\t"<<wordResult.pvalue<<"\t"<<FDR<<"\t"<<wordResult.expectedProb<<"\t"<<wordResult.expected<<"\t"<<wordResult.observedProb<<"\t"<<wordResult.observed<<"\t"<<outputWord<<endl;
						outputed++;
						cerr<<"t"<<endl;
					}
					


					if(FDRCutOff>=0 && FDR>FDRCutOff)
						break;
					
					if(outputed>=topCount)
					{
						cerr<<"top Count cut-offed:"<<(outputed)<<endl;
						break;
					}
					
				}
				
			
			
		}
	
	}

	cerr<<"u"<<endl;
	if(motifFam)
	{
		delete motifFam;
	}
	cerr<<"v"<<endl;

}




void printUsage(const char**argv)
{
	cerr<<argv[0]<<" -h filenamefg,filenamebg,headerRow1,startRow1,strandCol1,fromCol1,toCol1,N,CGBinInterval,FDRCutOff,pvalueCutOff,topCount,familyFile[.=no],noOverlap[=yes,no],usePSSM[=yes,no],PSSM.nShuffle,PSSM.FDR"<<endl;
}






int main(int argc,const char** argv)
{

	srand(time(NULL));

	//cerr<<StringUtil::replace("---AACT-CAG--","-","")<<endl;
	//return 1;
	/*string a="b       a";

	vector<string> splits;
	StringUtil::splitNoEmpty(a," ",splits);
	for(vector<string>::iterator i=splits.begin();i!=splits.end();i++)
		cerr<<(*i)<<endl;

	return 1;
	*/
	/*vector<int> input;
	vector<double> result;

	input.push_back(150);
	input.push_back(262);
	input.push_back(6543);
	input.push_back(18429);
	input.push_back(150);
	input.push_back(262);
	input.push_back(6543);
	input.push_back(-18429);
	JHypergeometric_Pvalue(input,result);

	for(vector<double>::iterator i = result.begin();i!=result.end();i++)
		cerr<<*i<<endl;


	return 1;*/

	/*MotifFamily motifFamily("./spliceAidFam.txt");
	PSSM* pssm=motifFamily.getPSSMOf("9G8");
	string seq("TAGGAACAGAGCGAGACGCAGACGAGAGAGACGAGATCTCCTCCTCGACGCGGAACTGCCTAAGAATGGGGGACGAGAGACGACGACGAGACGACGACGAGAGACGAGAGAGCACGACGCAGACTACTAGCTAGCATCAGCATCGACATCATCATCACAGCGACTACTACGACATCGACTAGCACTAGCATCACCTTACTATCTACGACATCTATCATCTATTTTTACGACATCACATCTACACCTTCTACTCTCTCCTCTCTCATCACTCATATCACTACGACTACGACTAC");
	cerr<<"seq length="<<seq.length()<<endl;
	pssm->setBackground(seq);
	cerr<<"bacground set"<<endl;
	//cerr<<"score at 1 "<<pssm->score(seq.c_str(),0)<<endl;

	pair<int,PSSM::ScoreSeqPosMap*> searchResult=pssm->getNumSigHitsFDR(seq,0.05,true,true,100);
	//int numSig=pssm->getNumSigHits(seq,0.05,100);
	int numSig=searchResult.first;
	PSSM::ScoreSeqPosMap *hits=searchResult.second;

	cerr<<"numSig="<<numSig<<endl;
	hits->printScorePosSeq(cerr,numSig);

	delete hits;

	if(pssm)
		delete pssm;

	return 1;*/

	cerr<<"MotifEnrichment Hypergeometric"<<endl;
	cerr<<__DATE__<<endl;
	if(argc<13)
		printUsage(argv);
	else
	{
		string method=argv[1];

		if(method=="-h")
		{
			if(argc<19)
				printUsage(argv);
			else
			{
				bool usePSSM=(string(argv[16])=="yes");
				int PSSMnShuffle=StringUtil::atoi(argv[17]);
				double PSSMFDRThreshold=StringUtil::atof(argv[18]);
				findEnrichedMotifByHGWE(argv[2],argv[3],StringUtil::atoi(argv[4]),StringUtil::atoi(argv[5]),StringUtil::atoi(argv[6]),StringUtil::atoi(argv[7]),StringUtil::atoi(argv[8]),StringUtil::atoi(argv[9]),StringUtil::atof(argv[10]),StringUtil::atof(argv[11]),StringUtil::atof(argv[12]),StringUtil::atoi(argv[13]),argv[14],(string(argv[15])=="yes"),usePSSM,PSSMnShuffle,PSSMFDRThreshold);
			}

		}
		else
		{
			cerr<<"unknown method "<<method<<endl;
		}
	}
	
	return 1;
}
