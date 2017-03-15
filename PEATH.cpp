#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>

#define EPSILON 0.00000000001

#define POPSIZE 100
#define OFFSIZE 50

using namespace std;

typedef struct {
	int fileNum;
	int fragNum;
}FileNum;

typedef struct {
	double D_hf;
	double D_hfstar;
}D_Entry;

typedef struct {
	char * str;
	double * qStr;
	double * qStr_star;
	int start_pos;
	int end_pos;
	int length;
}SubFragment;

typedef struct {
	SubFragment * subFragment;
	int subFragNum;
	int start_pos;
	int end_pos;
	int length;
}Fragment;

typedef struct {
	int start_frag;
	int end_frag;
	int problemSize;
}Block;

typedef struct {
	char * chr;
	double fitVal;
	double sum_D;
}Chromosome;

typedef struct {
	int start_frag;
	int end_frag;
	int endPos;
}FragsForPos;


int load_matrixFile();
double divideByBlock();
int haplotype_phasing(Block &pos, ofstream &);

Fragment * fragment;
FragsForPos * fragsForPos;

vector <Block> block;


int heteroNum;
int fragNum;
int realStartIndex;
char * tmpAns;

char matrixFileName[100];

													// uniform distribution function
unsigned seed = chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator(seed);
uniform_real_distribution<double> distribution(0.0, 1.0);



bool f(const FileNum&a, const FileNum&b) {
	return a.fragNum < b.fragNum;
}


bool compare_fit_val(Chromosome &left, Chromosome &right) {
	return left.fitVal > right.fitVal;
}

bool Greater2(Block &left, Block &right) {
	return left.problemSize > right.problemSize;
}


// hapNum에 해당하는 file의 frag data 가져오기
int load_matrixFile()
{
	FILE * matrixFile = fopen(matrixFileName, "r");

	fseek(matrixFile, 0, SEEK_END);
	int fileSize = ftell(matrixFile);
	fseek(matrixFile, 0, SEEK_SET);

	char buf[20];
	char * p;

	fgets(buf, 20, matrixFile);

	p=strtok(buf, " ");
	fragNum = atoi(p) -1;

	p=strtok('\0', " ");
	heteroNum = atoi(p);




	int fragDataSize = fileSize+fragNum+10;
	char * fragData = new char[fragDataSize];

	fread(fragData, sizeof(char), fragDataSize, matrixFile);	// 데이터 읽어오기

	fclose(matrixFile);

	fragment = new Fragment[fragNum];

	for (int i=0; i < fragNum; i++, p++)
	{
		if (i==0)
			p = strtok(fragData, " ");	// subFragment 갯수
		else
			p = strtok(p, " ");	// subFragment 갯수

		int subFragNum = atoi(p); // i번째 Fragment의 subFragment 갯수

		fragment[i].subFragNum= subFragNum;
		fragment[i].subFragment = new SubFragment[subFragNum]; // i번쨰 fragment의 Subfragment갯수만큼 할당

		p = strtok('\0', " "); // 이름 버림

							   // 데이터 받아오는부분 -1로해서 받아오기

		for (int j=0; j < subFragNum; j++)
		{
			p = strtok('\0', " ");
			fragment[i].subFragment[j].start_pos = atoi(p)-1; // start 갱신

			p = strtok('\0', " ");
			fragment[i].subFragment[j].str = p;
			fragment[i].subFragment[j].length = strlen(fragment[i].subFragment[j].str);
			fragment[i].subFragment[j].end_pos = fragment[i].subFragment[j].start_pos + fragment[i].subFragment[j].length -1;// 길이와 end 갱신


		} // end j

		fragment[i].start_pos = fragment[i].subFragment[0].start_pos;
		fragment[i].end_pos = fragment[i].subFragment[subFragNum-1].end_pos;
		fragment[i].length = fragment[i].end_pos - fragment[i].start_pos +1; // length 갱신

		p = strtok('\0', "\n"); // " " -> "\n" 으로 수정

		double * qStr = new double[fragment[i].length];
		double * qStr_star = new double[fragment[i].length];

		for (int j=0; j<subFragNum; j++)
		{
			for (int k=0; k<fragment[i].subFragment[j].length; k++, p++) {
				int qual_j = (int)(*p);
				qual_j -=33;
				qStr[k] = pow(10.0, ((-1*(double)qual_j)/10.0));
				qStr_star[k] = 1-qStr[k];
			}

			fragment[i].subFragment[j].qStr = qStr;
			fragment[i].subFragment[j].qStr_star = qStr_star;

			qStr += (fragment[i].subFragment[j].length);	// 여기서 p가 바뀜
			qStr_star += (fragment[i].subFragment[j].length);	// 여기서 p가 바뀜

		} // end for j

	} // end i


	  // fragment 잘 받아왔는지 파일로 확인/////////////////////////////
	ofstream fragmentOutFile;
	fragmentOutFile.open("fragmentOut.txt");

	for (int i=0; i<fragNum; i++) {
		fragmentOutFile<< "start: " << fragment[i].start_pos << " ";
		fragmentOutFile<< "end: " << fragment[i].end_pos << " ";
		fragmentOutFile<< "length: " << fragment[i].length << " ";
		fragmentOutFile<< "subFragNum: " << fragment[i].subFragNum <<endl;

		for (int j=0; j<fragment[i].subFragNum; j++) {
			fragmentOutFile<< "start: " << fragment[i].subFragment[j].start_pos << " end: " << fragment[i].subFragment[j].end_pos << " length: " << fragment[i].subFragment[j].length << " str:" << fragment[i].subFragment[j].str << " qual:" << fragment[i].subFragment[j].qStr << endl;
		}
		fragmentOutFile << "==============================================================" << endl;
	}

	fragmentOutFile.close();
	//////////////////////////////////////////////////////////////////


	return 0;
}

// 해당 block의 phasing 작업
void calc_fitVal(Chromosome& chromosome, Block& block)
{
	//call_count_GA++;

	int startFrag, endFrag;
	double sum_D=0.0;

	startFrag = block.start_frag;
	endFrag=block.end_frag;

	realStartIndex = fragment[block.start_frag].start_pos;

	for (int i=startFrag; i <= endFrag; i++) {

		double D_hf=0.0, D_hstarf=0.0;

		for (int j=0; j<fragment[i].subFragNum; j++) {

			for (int k=0; k<fragment[i].subFragment[j].length; k++) {

				double q_j = fragment[i].subFragment[j].qStr[k];
				double q_j_star = fragment[i].subFragment[j].qStr_star[k];

				if (chromosome.chr[fragment[i].subFragment[j].start_pos+k-realStartIndex] != fragment[i].subFragment[j].str[k]) {
					D_hf += q_j_star;
					D_hstarf += q_j;
				}
				else {
					D_hf += q_j;
					D_hstarf += q_j_star;
				}

			}// end for k (작은 서브프래그먼트가 끝난시점)

		} // j

		if (D_hf < D_hstarf)
			sum_D += D_hf;
		else
			sum_D += D_hstarf;

	}// i

	chromosome.sum_D=sum_D;
	chromosome.fitVal=1.0/sum_D*1000000;
}


void calc_fitVal_rangeSwitch(D_Entry *& D, Chromosome& chromosome, Block& block, int pos)
{
	//call_count_range++;

	int startFrag, endFrag;
	double sum_D=0.0;

	//if (fragsForPos[pos+fragment[block.start_frag].start_pos].start_frag == -2) return;// population.fitVal;
	if (fragsForPos[pos+fragment[block.start_frag].start_pos].start_frag > fragsForPos[pos+fragment[block.start_frag].start_pos].end_frag) return;// population.fitVal;

	sum_D=chromosome.sum_D;

	pos = pos+fragment[block.start_frag].start_pos;
	startFrag = fragsForPos[pos].start_frag;
	endFrag = fragsForPos[pos].end_frag;

	realStartIndex = fragment[block.start_frag].start_pos;

	// mod for fragments including pos
	for (int i=startFrag; i <= endFrag; i++) {

		double D_hf=D[i].D_hf, D_hfstar=D[i].D_hfstar;

		if (D_hf < D_hfstar) sum_D -= D_hf;
		else sum_D -= D_hfstar;

		for (int j=0; j<fragment[i].subFragNum; j++) {

			if (pos > fragment[i].subFragment[j].end_pos)
				continue;
			if (pos < fragment[i].subFragment[j].start_pos)
				break;

			int k=pos-fragment[i].subFragment[j].start_pos;

			double q_j = fragment[i].subFragment[j].qStr[k];
			double q_j_star = fragment[i].subFragment[j].qStr_star[k];

			// mod 설명, chr의 pos에 해당하는 비트가 flip됬다고 가정하고 반대로 계산, != -> ==			
			if (chromosome.chr[pos-fragment[block.start_frag].start_pos] == fragment[i].subFragment[j].str[k]) {
				D_hf += q_j_star-q_j;
				D_hfstar+= q_j-q_j_star;
			}
			else {
				D_hf += q_j - q_j_star;
				D_hfstar+= q_j_star-q_j;
			}


		} // end for j : 모든 subFragment가 끝난 시점




		  // 모든 i마다 수행됨
		if (D_hf < D_hfstar)
			sum_D += D_hf;
		else
			sum_D += D_hfstar;

		D[i].D_hf = D_hf;
		D[i].D_hfstar = D_hfstar;

	}// end for i : 모든 Fragment 끝난 시점

	chromosome.sum_D=sum_D;
	chromosome.fitVal=1.0/sum_D*1000000;
}


double calc_fitVal_singleSwitch(D_Entry *& D, Chromosome& chromosome, Block& block, int pos)
{
	//call_count_single++;

	int startFrag, endFrag;
	double sum_D=0.0;

	//if (fragsForPos[pos+fragment[block.start_frag].start_pos].start_frag == -2) return chromosome.sum_D;// population.fitVal;
	if (fragsForPos[pos+fragment[block.start_frag].start_pos].start_frag > fragsForPos[pos+fragment[block.start_frag].start_pos].end_frag) return chromosome.sum_D;// population.fitVal;

	sum_D=chromosome.sum_D;

	pos = pos+fragment[block.start_frag].start_pos;
	startFrag = fragsForPos[pos].start_frag;
	endFrag = fragsForPos[pos].end_frag;

	realStartIndex = fragment[block.start_frag].start_pos;

	// mod for fragments including pos
	for (int i=startFrag; i <= endFrag; i++) {

		double D_hf=D[i].D_hf, D_hfstar=D[i].D_hfstar;

		if (D_hf < D_hfstar) sum_D -= D_hf;
		else sum_D -= D_hfstar;

		for (int j=0; j<fragment[i].subFragNum; j++) {

			if (pos > fragment[i].subFragment[j].end_pos)
				continue;
			if (pos < fragment[i].subFragment[j].start_pos)
				break;

			int k=pos-fragment[i].subFragment[j].start_pos;

			double q_j = fragment[i].subFragment[j].qStr[k];
			double q_j_star = fragment[i].subFragment[j].qStr_star[k];

			// mod 설명, chr의 pos에 해당하는 비트가 flip됬다고 가정하고 반대로 계산, != -> ==
			if (chromosome.chr[pos-fragment[block.start_frag].start_pos] == fragment[i].subFragment[j].str[k]) {
				D_hf += q_j_star-q_j;
				D_hfstar+= q_j-q_j_star;
			}
			else {
				D_hf += q_j - q_j_star;
				D_hfstar+= q_j_star-q_j;
			}

		} // end for j : 모든 subFragment가 끝난 시점


		  // 모든 i마다 수행됨
		if (D_hf < D_hfstar)
			sum_D += D_hf;
		else
			sum_D += D_hfstar;


	}// end for i : 모든 Fragment 끝난 시점

	return sum_D;
}

// mod 이름 변경
double calc_fitVal_D(D_Entry *& D, Chromosome& chromosome, Block& block)
{
	int startFrag, endFrag;
	double sum_D=0.0;

	startFrag = block.start_frag;
	endFrag=block.end_frag;

	realStartIndex = fragment[block.start_frag].start_pos;

	// mod for all fragments in block
	for (int i=startFrag; i <= endFrag; i++) {

		double D_hf=0.0, D_hstarf=0.0;

		for (int j=0; j<fragment[i].subFragNum; j++) {

			for (int k=0; k<fragment[i].subFragment[j].length; k++) {


				double q_j = fragment[i].subFragment[j].qStr[k];
				double q_j_star = fragment[i].subFragment[j].qStr_star[k];

				if (chromosome.chr[fragment[i].subFragment[j].start_pos+k-realStartIndex] != fragment[i].subFragment[j].str[k]) {
					D_hf += q_j_star;
					D_hstarf += q_j;
				}
				else {
					D_hf += q_j;
					D_hstarf += q_j_star;
				}

			}// end for k (작은 서브프래그먼트가 끝난시점)


		} // j

		if (D_hf < D_hstarf)
			sum_D += D_hf;
		else
			sum_D += D_hstarf;

		D[i].D_hf = D_hf;
		D[i].D_hfstar = D_hstarf;

	}// i

	chromosome.sum_D=sum_D;
	chromosome.fitVal=1.0/sum_D*1000000;
	return sum_D;
}


// hapNum에 해당하는 file을 phasing하는 작업
int procedure() {

	load_matrixFile(); // fragNum, heteroNum 갱신

					   // fragment 잘 받아왔는지 파일로 확인/////////////////////////////
	ofstream fragmentOutFile;
	fragmentOutFile.open("fragmentOut.txt");

	for (int i=0; i<fragNum; i++) {
		fragmentOutFile<< "start: " << fragment[i].start_pos << " ";
		fragmentOutFile<< "end: " << fragment[i].end_pos << " ";
		fragmentOutFile<< "length: " << fragment[i].length << " ";
		fragmentOutFile<< "subFragNum: " << fragment[i].subFragNum <<endl;

		for (int j=0; j<fragment[i].subFragNum; j++) {
			fragmentOutFile<< "start: " << fragment[i].subFragment[j].start_pos << " end: " << fragment[i].subFragment[j].end_pos << " length: " << fragment[i].subFragment[j].length << " str:" << fragment[i].subFragment[j].str << " qual:" << fragment[i].subFragment[j].qStr;
		}
		fragmentOutFile << "==============================================================" << endl;
	}

	fragmentOutFile.close();
	// fragment 잘 받아왔는지 파일로 확인////////////////////////////////////

	fragsForPos = new FragsForPos[heteroNum]; // 0으로 초기화
	for (int i=0; i<heteroNum; i++)
		fragsForPos[i]={ -1,-1 };


	//////////////////////////////////////////////////////////

	// 초기화 후 block 별로 나누는 작업
	divideByBlock();

	printf("block# :%d\n", block.size());

	ofstream phasingResultFile;

	char phasingResFileName[101];

	strcpy(phasingResFileName, matrixFileName);
	strcat(phasingResFileName, "_Result.txt"); // 여기까지 수정

	phasingResultFile.open(phasingResFileName);
	phasingResultFile.precision(5);

	int totalBlockSize=0;

	// block 별로 Phasing

	for (int i=0; i<block.size(); i++) {
			
		printf("Block # %d\n", i+1);
		printf("\n Block Range : %d %d\n", block[i].start_frag, block[i].end_frag);
		printf("\n Block Length : %d\n", block[i].problemSize);
			
		phasingResultFile << "=========================================================\nBlock Number : " << i+1 << '\n';

		totalBlockSize+=block[i].problemSize;

		haplotype_phasing(block[i], phasingResultFile);
	}

	phasingResultFile << "Total Block Problem Size : " << totalBlockSize << '\n';

	phasingResultFile.close();

	delete fragment;

	return 0;
}


void GA(Block &block, Chromosome & genAns) {

	// mod cnt 배열생성

	bool isStop = true;
	double maxFit=-1;
	int stopCnt=0;

	//int fragmentNumber = pos.end-pos.start+1;
	Chromosome * population = new Chromosome[POPSIZE];
	Chromosome * offspring = new Chromosome[OFFSIZE];

	int * cnt = new int[block.problemSize]{ 0 };
	int * cnt2 = new int[block.problemSize]{ 0 };

	// population assignment, fitness value도 설정해야함
	for (int i=0; i<POPSIZE; i++) {
		population[i].chr = new char[block.problemSize+1];
		population[i].fitVal = 0;
		population[i].sum_D=0;
		population[i].chr[block.problemSize] = '\0';
	}

	// mod offspring 따로 생성 X
	// offspring assignment
	for (int i=0; i<OFFSIZE; i++) {
		offspring[i].chr = new char[block.problemSize+1];
		offspring[i].fitVal = 0;
		offspring[i].sum_D=0;
		offspring[i].chr[block.problemSize]= '\0';
	}

	// 랜덤으로 population init
	for (int i=0; i<POPSIZE; i++) {
		for (int j=0; j<block.problemSize; j++) {

			if (distribution(generator)<0.5)
				population[i].chr[j]='0';
			else {
				population[i].chr[j]='1';

				if (i<OFFSIZE)
					cnt[j]++;
				else
					cnt2[j]++;
			}
		}
	}

	// 초기해집단 적합도 평가, mod 위에서 init하면서 같이
	for (int i=0; i<POPSIZE; i++) {
		calc_fitVal(population[i], block);
	}

	sort(population, &population[POPSIZE], compare_fit_val);


	// GA simulation
	int iter = 50;
	while (iter--) {	// mod 수렴하는 정도에 따라서 조절

		double prob;
		for (int i=0; i<block.problemSize; i++)
			cnt[i]+=cnt2[i];

		for (int i = OFFSIZE; i<POPSIZE-1; i++) {
			for (int j=0; j<block.problemSize; j++) {

				prob = (double)cnt[j]/POPSIZE;

				if (population[i].chr[j] == '1')
					cnt2[j]--;

				if (distribution(generator)<prob) {
					population[i].chr[j] = '1';
					cnt2[j]++;
				}
				else
					population[i].chr[j] = '0';

			}

		}


		for (int i=OFFSIZE; i<POPSIZE-1; i++)
			calc_fitVal(population[i], block);

		// 정렬
		sort(population, &population[POPSIZE], compare_fit_val);


		if (population[0].fitVal - maxFit > EPSILON) {
			maxFit = population[0].fitVal;
			stopCnt = 0;
		}
		else
			stopCnt++;

		if (stopCnt >= 10)
			break;
	}

	genAns.chr = new char[block.problemSize+1];
	strcpy(genAns.chr, population[0].chr);
	genAns.fitVal = population[0].fitVal;
	genAns.sum_D = population[0].sum_D;

}


void range_switch_procedure(D_Entry *& D, Block &block, Chromosome & genAns) {

	double maxSum = genAns.sum_D;

	char * tmpStr = new char[block.problemSize+1];
	tmpStr[block.problemSize] = '\0';
	strcpy(tmpStr, genAns.chr);

	int maxIdx = -1;
	bool isImp = true;

	calc_fitVal_D(D, genAns, block); // mod

	while (isImp) {

		isImp = false;

		for (int i=0; i<block.problemSize; i++) {

			genAns.chr[i] == '0' ? (genAns.chr[i] = '1') : (genAns.chr[i] = '0');
			calc_fitVal_D(D, genAns, block);


			if ((maxSum - genAns.sum_D) > EPSILON) {

				isImp = true;
				maxSum = genAns.sum_D;
				maxIdx = i;
			}

		}


		// 0~maxIdx까지 뒤집기
		if (isImp) {
			for (int i=0; i<=maxIdx; i++)
				genAns.chr[i] == '0' ? (genAns.chr[i] = '1') : (genAns.chr[i] = '0');

			for (int i=0; i<=maxIdx; i++)
				tmpStr[i] == '0' ? (tmpStr[i] = '1') : (tmpStr[i] = '0');

			calc_fitVal_D(D, genAns, block); // D랑은 상관이 없음 mod calc_fitval_D mod
			genAns.sum_D = maxSum;
		}

	}
	delete tmpStr;
}


void single_switch_procedure(D_Entry *& D, Block &block, Chromosome & genAns) {

	// 임의로 추가

	char * tmpStr = new char[block.problemSize+1];
	tmpStr[block.problemSize] = '\0';
	strcpy(tmpStr, genAns.chr);

	double maxSum = genAns.sum_D;

	int maxIdx = -1;
	bool isImp = true;

	isImp = true;
	while (isImp) {

		double tmpSum_D;
		isImp = false;

		for (int i=0; i<block.problemSize; i++) {

			tmpSum_D = calc_fitVal_singleSwitch(D, genAns, block, i);	// 평가

			if (maxSum > tmpSum_D) {
				isImp = true;
				maxIdx = i;
				maxSum = tmpSum_D;
			}

		}

		if (isImp) {
			genAns.chr[maxIdx] == '0' ? (genAns.chr[maxIdx] = '1') : (genAns.chr[maxIdx] = '0');	// 최종적으로 한비트 바꿈
			calc_fitVal_D(D, genAns, block);

		}

	}
}

int haplotype_phasing(Block &block, ofstream & phasingResultFile) {

	Chromosome genAns;
	GA(block, genAns);


	D_Entry * D = new D_Entry[fragNum];
	calc_fitVal_D(D, genAns, block);


	range_switch_procedure(D, block, genAns);

	single_switch_procedure(D, block, genAns);

	range_switch_procedure(D, block, genAns);
	single_switch_procedure(D, block, genAns);



	// 최종적으로 가장 위에있는 (가장 좋은) 해 출력


	printf("%s\n", genAns.chr);
	phasingResultFile << genAns.chr << endl;

	phasingResultFile << '\n';
	phasingResultFile << "Block Length : " << block.problemSize << "\n";


	delete D;
	return 0;
}


double divideByBlock()
{
	int cutStartFrag=0;	// 자르기 시작하는 fragment 번호

	int maxEndIdx = fragment[0].end_pos;
	bool chk = false;
	for (int i=1; i<fragNum; i++) {	// 2번째부터 시작

		if (fragment[i].start_pos > maxEndIdx) {
			block.push_back({ cutStartFrag,i-1, maxEndIdx - fragment[cutStartFrag].start_pos + 1 }); // 구간 저장하고
			cutStartFrag = i; // 새로운 컷 번호 갱신

		}
		maxEndIdx = max(maxEndIdx, fragment[i].end_pos); // 최대 인덱스 유지
	}

	// 마지막 무조건 남게 되어있음
	block.push_back({ cutStartFrag, fragNum-1, maxEndIdx - fragment[cutStartFrag].start_pos + 1 }); // 구간 저장하고

	ofstream blockFile;
	blockFile.open("BlockPos.txt");
	for (int i=0; i<block.size(); i++) {
		blockFile << "[" << i+1 << "] " << block[i].start_frag << " " << block[i].end_frag << " " << block[i].problemSize << endl;

	}
	blockFile.close();


	// start_frag갱신

	int cutStart=0;
	for (int i=0; i<fragNum; i++) {

		for (; cutStart <= fragment[i].end_pos; cutStart++)
			if (fragsForPos[cutStart].start_frag == -1)
				fragsForPos[cutStart].start_frag = i;

	}

	// end_frag갱신

	int cutEnd = heteroNum-1;
	for (int i=fragNum-1; i>=0; i--) {

		for (; cutEnd >= fragment[i].start_pos; cutEnd--)
			if (fragsForPos[cutEnd].end_frag == -1)
				fragsForPos[cutEnd].end_frag = i;

	}


	ofstream rangePosFile;
	rangePosFile.open("rangePos.txt");
	for (int i=0; i<heteroNum; i++) {
		rangePosFile << "[" << i << "] " << fragsForPos[i].start_frag << " " << fragsForPos[i].end_frag <<endl;
	}
	rangePosFile.close();



	ofstream fragmentForm("FragmentForm.SORTED");

	for (int i=0; i<block.size(); i++) {

		fragmentForm << "block # : " << i+1 << endl;

		for (int j= block[i].start_frag; j<=block[i].end_frag; j++) {

			int startIdx = fragment[block[i].start_frag].start_pos;

			char * frag  = new char[block[i].problemSize+1];
			for (int k=0; k<block[i].problemSize; k++)
				frag[k] = '-';
			frag[block[i].problemSize] = '\0';


			for (int k = 0; k<fragment[j].subFragNum; k++) {

				int idx = fragment[j].subFragment[k].start_pos - startIdx;

				for (int l=0; l<fragment[j].subFragment[k].length; l++) {


					frag[idx+l] = fragment[j].subFragment[k].str[l];
				}

			}
			fragmentForm << frag << endl;
		}

	}

	fragmentForm.close();

	return 1;
}


// 프로그램 input으로 fragment file이름 전달 ex) chr22_matrix
int main(int argc, char ** argv) {

	printf("Enter the Matrix File Name: ");
	scanf("%s", matrixFileName);
	int chrome_num=0;

	procedure();

	// matrix file 없을시 에러 처리

	return 0;
}

