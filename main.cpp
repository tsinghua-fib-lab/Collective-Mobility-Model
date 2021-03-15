#include <iostream>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include <time.h>
#include <algorithm>
#include <map>
#include <stdlib.h>
#include "drand48.h"
#include <ppl.h>
#include <string>
#include <windows.h>
#include <hash_map>
#include <bitset>
#include <direct.h>

#pragma hdrstop

using namespace std;
using namespace Concurrency;

typedef unsigned int uint;
#define const_a 0x5DEECE66DLL
#define const_c 0xB
#define DRAND48_MAX 0xFFFFFFFFFFFFLL
#define LRAND_MAX 0xFFFFFFFFL

#ifdef __BORLANDC__
inline double drand48()
{
	return double(_lrand())/LRAND_MAX;
}
#endif

const int L = 300, A = L*L;//***simulation space size
int N;
typedef int myvec;

int p1[A], p2[A], p3[A], p4[A];

int *pop = p1, *pre_pop = p2, *avg_pop = p3;
int *pop_static=p4;

combinable<int> pop_temp[A];
combinable<double> pa_time;
combinable<double> ex_time;
combinable<int> pa_num;
combinable<int> ex_num;

//model parameters

double delta=1;
double rho0=0; 
double alpha=-2.55;
double gamma= 0.21;
int ** snap_shots;

//simulation system parameters
int memory_flag = 0;
unsigned int backet_num=5000;
int pop_threshold = 30;
int window_output = 50;
int window_history = 1;
const int pop_sample_num = 100;

const int dist_threshold = L/2;
const int candidate_size = 4*dist_threshold*dist_threshold;

class AliasTable
{
	public:

		void deletTable()
		{
		free(alias);
		free(label);
		free(prob);				
		}

		void InitAliasTable(int prob_length, double *desire_prob, int *tag)
		{
			block_num = prob_length;
			alias = (int *)malloc(block_num*sizeof(int));
			label = (int *)malloc(block_num*sizeof(int));
			prob = (double *)malloc(block_num*sizeof(double));
			if (alias == NULL || prob == NULL || label == NULL)
			{
				printf("Error: memory allocation failed!\n");
				exit(1);
			}

			double *norm_prob = (double*)malloc(block_num*sizeof(double));
			int *large_block = (int*)malloc(block_num*sizeof(int));
			int *small_block = (int*)malloc(block_num*sizeof(int));
			if (norm_prob == NULL || large_block == NULL || small_block == NULL)
			{
				printf("Error: memory allocation failed!\n");
				exit(1);
			}

			double sum = 0;
			int cur_small_block, cur_large_block;
			int num_small_block = 0, num_large_block = 0;

			for (int k = 0; k != block_num; k++) sum += desire_prob[k];//ͳ�����бߵ�Ȩ�غ�
			for (int k = 0; k != block_num; k++) norm_prob[k] = desire_prob[k] * block_num / sum;//ÿ����ռƽ��Ȩ�ض���

			for (int k = block_num - 1; k >= 0; k--)//��׼��Ȩ�ط�Ϊ����
			{
				if (norm_prob[k]<1)
					small_block[num_small_block++] = k;
				else
					large_block[num_large_block++] = k;
			}

			int k = 0;
			while (num_small_block && num_large_block)
			{
				cur_small_block = small_block[--num_small_block];  //cur_small_block�ǵ�ǰС�����
				cur_large_block = large_block[--num_large_block];  //cur_large_block�ǵ�ǰ��ߵ����

				prob[k] = norm_prob[cur_small_block];//����һ�����ʸ���prob����
				label[k] = tag[cur_small_block];
				alias[k] = tag[cur_large_block];          //��ǰС�߶�Ӧ�ı�������ֵ�ǵ�ǰ��ߵ����
				norm_prob[cur_large_block] = norm_prob[cur_large_block] + norm_prob[cur_small_block] - 1;

				if (norm_prob[cur_large_block] < 1)
					small_block[num_small_block++] = cur_large_block;
				else
					large_block[num_large_block++] = cur_large_block;
				k++;
			}

			while (num_large_block) 
			{
				prob[k] = 1;
				alias[k] = -1;
				label[k] = tag[large_block[--num_large_block]]; 
				k++;
			}
			while (num_small_block) 
			{
				prob[k] = 1;
				alias[k] = -1;
				label[k] = tag[small_block[--num_small_block]];
				k++;
			}

			free(norm_prob);
			free(small_block);
			free(large_block);
		}

		int GenerateSamples(double rand_value1, double rand_value2)
		{
			int k = (int)block_num * rand_value1;
			return rand_value2 < prob[k] ? label[k] : alias[k];
		}

		void PrintTable()
		{
			ofstream file("test.dat");
		    for (int i = 0; i != block_num; i++) file<<"("<<prob[i]<<","<<alias[i]<<")";//printf("%f ", prob[i]);
				file<<endl;//printf("\n");   
			//for (int i = 0; i != block_num; i++) file<<alias[i]<<" ";//printf("%d ", alias[i]);
				//file<<endl;//printf("\n");		
		}

	private:
		int *alias;
		int *label;
		double *prob;
		int block_num;
};

map<myvec, AliasTable> tables;

class Regions
{
public:
   double dist[candidate_size];
   int  x[candidate_size];
   int  y[candidate_size];
   int  cur;
};

Regions cutoffRegion;

class SortedArray
{
	private:
		int *label;
		int *frequency;
		int capacity;
		int cur;

	public:
		SortedArray(int len)
		{
			this->capacity = len;
			this->frequency = (int *) malloc(len*sizeof(int));
			this->label = (int *) malloc(len*sizeof(int));
			this->cur = 0;
			if (this->frequency==NULL||this->label==NULL)
				cout<<"Fail to allocate memory!"<<endl;

		}

		SortedArray()
		{
		}

		~SortedArray()
		{
			free(this->frequency);
			free(this->label);
		}

		void insert(int element_label)
		{
			if(cur<capacity)
			{
				this->label[cur] = element_label;
				this->frequency[cur] = 1;
				cur++;
			}
			else
			{
				cout<<"Insert error. Out of memory!"<<endl;
			}
		}

		void select(int index, int &output_label, int &output_frequency)
		{
			if (index>this->cur)
			{
				cout<<"Select out of index!"<<endl;
			}
			else
			{
				output_label = this->label[index];
				output_frequency = this->frequency[index];
			}
		}

		void add(int index)
		{
			if (index>this->cur)
			{
				cout<<"Add out of index!"<<endl;
			}
			else
			{
				if (index==0)
				{
					this->frequency[0]++;
				}
				else
				{
					int label_temp = this->label[index];
					int freq_temp = this->frequency[index]+1;
					int index_temp=index;
					while((freq_temp>this->frequency[index_temp-1])&&index_temp>0)
					{
						this->frequency[index_temp] = this->frequency[index_temp-1];
						this->label[index_temp] = this->label[index_temp-1];
						index_temp--;
					}

					this->frequency[index_temp] =  freq_temp;
					this->label[index_temp] = label_temp;
				}				
			}
		}

		void PrintTable()
		{
			//ofstream file("test.dat");
			for (int i = 0; i != capacity; i++) cout<<"("<<this->frequency[i]<<","<<this->label[i]<<")";//printf("%f ", prob[i]);
				cout<<endl;
		}
};

inline myvec vec2(int x, int y) 
{
	return x + y*L;
}

void initCutoffRegion()
{
	int disttemp = 0;
	int count = 0;
	for(int tempx = - dist_threshold; tempx <= dist_threshold; tempx++)
	{
		for(int tempy = - dist_threshold; tempy <= dist_threshold; tempy++)
		{
			if(tempx==0 && tempy==0)
				continue;

			disttemp = sqrt((double)tempx*tempx+tempy*tempy);

			if(disttemp>dist_threshold)
				continue;

			cutoffRegion.dist[count] = disttemp;
			cutoffRegion.x[count] = tempx;
			cutoffRegion.y[count] = tempy;
			count ++;
			cutoffRegion.cur = count;
		}	
	}

}

class t_user
{
	bitset <A> unvisited_flag;
	int visited_num;
	myvec pos;
	long long seed_pre;
	long times;
	SortedArray * loc_freq;

	void explore(int t)/* explore new urban */
	{
		myvec npos;/* next move */
		int count_flag = 1;
		if (tables.find(this->pos)!=tables.end())
		{
			int count = 0;
			while(count<100)
			{   count ++;
				npos=tables[pos].GenerateSamples(this->drand48(),this->drand48());
				if(this->unvisited_flag[npos])
					break;
			}
			if (count<100)
			{
				count_flag = 0;
			}
		}

		if(count_flag)
		{
			double *prob;
			int *postemp;
			prob = (double *) malloc(candidate_size*sizeof(double));
			postemp = (int *) malloc(candidate_size*sizeof(int));
			double sum = 0;
			int tempx = pos%L;
			int tempy = pos/L;
			int x = 0;
			int y = 0;
			int count = 0;

			for(int i = 0;i<cutoffRegion.cur;i++)
			{
				x = tempx + cutoffRegion.x[i];
				y = tempy + cutoffRegion.y[i];

				if(x>=L)
					x -= L;
				else if(x<0)
					x += L;

				if(y>=L)
					y -= L;
				else if(y<0)
					y += L;

				if(vec2(x,y) == pos)
					continue;
				
				if(!this->unvisited_flag[vec2(x,y)])
					continue;
				prob[count] = (pre_pop[vec2(x, y)]+rho0)*pow(cutoffRegion.dist[i], alpha);
				sum += prob[count];
				postemp[count] = vec2(x, y);
				count ++;
			}

			if(sum==0)
			{
				count = 0;
				for(int i = 0;i<cutoffRegion.cur;i++)
				{
					x = tempx + cutoffRegion.x[i];
					y = tempy + cutoffRegion.y[i];

					if(x>=L)
						x -= L;
					else if(x<0)
						x += L;

					if(y>=L)
						y -= L;
					else if(y<0)
						y += L;

					if(vec2(x,y) == pos)
						continue;

					if(!this->unvisited_flag[vec2(x,y)])
						continue;
					prob[count] = pow(cutoffRegion.dist[i], alpha);
					sum += prob[count];
					postemp[count] = vec2(x, y);
					count ++;
				}				
			}
			double cut = this->drand48()*sum;
			sum = 0;
			for(int i = 0;i<count;i++)
			{
				sum += prob[i];
				if(sum>cut)
				{
					npos = postemp[i];
					break;
				} 
			}
			free(prob);
			free(postemp);
		}
    if(memory_flag)
	{
		this->unvisited_flag[npos]=false; 
		this->loc_freq->insert(npos);
		this->visited_num++;
		set_pos(npos);/* people move from pos to npos,change the population array */
		this->times++;
	}
	else
	{
		set_pos(npos);/* people move from pos to npos,change the population array */
		this->times++;
	}
	}

	void pa_return()/* move/return to an old location */
	{
		int  sum = 0;
		int  freq = 0;
		int  npos = -1;
		int  index;
		int cut = int(this->drand48()*this->times);

		for (index = 0; index<this->visited_num;  index++)
		{
			this->loc_freq->select(index, npos, freq);
			sum += freq;
			if (sum >= cut) break;/* >threshold */
		}

		this->loc_freq->add(index);

		set_pos(npos);
		this->times++;
	}

	void set_pos(const myvec& x) 
	{
		pop_temp[pos].local() --;
		this->pos = x;
		pop_temp[pos].local() ++;
	}
	
public:
	t_user() 
	{

	}

	void init(const myvec& p)/* p is the source point */
	{
		this->pos = p;	
		this->visited_num = 1;
		for (int i=0; i< A; i++)
			this->unvisited_flag[i]=true; 
		if(memory_flag) this->unvisited_flag[p]=false; 
		this->times = 1;
		this->loc_freq = new SortedArray(backet_num);
		this->loc_freq->insert(p);
		pop[p] ++;		
	}

	void run_one_step(int t)
	{
		clock_t temp;
		if ( this->drand48() < delta*pow(this->visited_num, -gamma))/* the case of exploring new urban,probability is delta*pow(N, -gamma)) */
			{
				temp = clock();
				explore(t);
				ex_time.local()  += (double)(clock()-temp)/CLOCKS_PER_SEC;
				ex_num.local()++;
			}
		else 
			{
				temp = clock();
				pa_return();
				pa_time.local()  += (double)(clock()-temp)/CLOCKS_PER_SEC;
				pa_num.local()++;
			}
	}

	void srand48(long int seedval)
	{
		seed_pre = ((long long int)seedval)<<16 | 0x00000000330ELL;
	}

	double drand48()
	{
		seed_pre = (const_a * seed_pre + const_c) % (DRAND48_MAX+1);
		unsigned long int seed = seed_pre >> 16;
		return seed / (LRAND_MAX + 1.);
	}

	unsigned long int lrand48()
	{
		seed_pre = (const_a * seed_pre + const_c) % (DRAND48_MAX+1);
		return seed_pre >> 16;
	}

};

class t_sta
{
	int count[pop_sample_num];
	int area[pop_sample_num];
	int x;
	int y;

	public:
	t_sta() 
	{ 
		fill_n(count, pop_sample_num,0); 
		fill_n(area, pop_sample_num,0); 
	}

	void init(int x, int y)
	{
		fill_n(count, pop_sample_num,0); 
		fill_n(area, pop_sample_num,0); 
		this->x = x;
		this->y = y;
	}

	void add(int target_x, int target_y, int p)
	{
		double r = sqrt((double)(target_x-x)*(target_x-x)+(target_y-y)*(target_y-y));
		int i = pop_sample_num*r/(L/2); 
		if (i >= pop_sample_num) return;
		count[i] += p; /* total population where dist = r */
		area[i] += 1;
	}

	void save(const char* name1)
	{
		ofstream s(name1,ios::app);

		for (int i = 0; i < pop_sample_num; i ++)
		{
			if(area[i]==0||count[i]==0)
			{
				s << "0\t";
			}
			else
			{
				s << float(count[i])/area[i] << '\t';
			}
		}
		s<<endl;
		s.close();
	}
};


struct t_neighbour : public vector<myvec>	
{
	int flag;	
	int pop_unit;  

	t_neighbour() : flag(-1),pop_unit(0){ } 

	void set_flag(map<myvec,t_neighbour>& u, int f)
	{
		flag = f;
		vector<myvec> candidate = *this;
		while(!candidate.empty())
		{
			myvec p = candidate.front();/* stack--LIFO */
			candidate.front() = candidate.back();
			candidate.pop_back();

			t_neighbour& n = u[p];
			if(n.flag == f)
				continue;
			n.flag = f;

			for (int i = 0; i < n.size(); i ++)
			{
				myvec p1 = n[i];/* n[i].first */
				if (u[p1].flag == -1) candidate.push_back(p1);
			}
		}
	}

};

typedef map<myvec,t_neighbour> urban_type;


void urban_size(const char* name, const char* name2, urban_type& u)
{
	for (urban_type::iterator i = u.begin(); i != u.end(); ++ i)
	{
		int x = i->first%L,  y = i->first/L, c;
		c = vec2((x+1)%L, y); if (pop[c]) i->second.push_back(c);
		c = vec2((x+L-1)%L, y); if (pop[c]) i->second.push_back(c);
		c = vec2(x, (y+1)%L); if (pop[c]) i->second.push_back(c);
		c = vec2(x, (y+L-1)%L); if (pop[c]) i->second.push_back(c);
	}

	int flag = 0;
	for (urban_type::iterator i = u.begin(); i != u.end(); ++ i)
	{
		if (i->second.flag != -1) continue;
		i->second.set_flag(u, flag);
		flag ++;
	}

	vector<int> size(flag, 0);
	vector<int> u_pop(flag,0);

	for (urban_type::iterator i = u.begin(); i != u.end(); ++ i)
	{
		int x = i->first%L;
		int y =  i->first/L;
		size[i->second.flag] ++;
		u_pop[i->second.flag] += i->second.pop_unit;
	}
	ofstream ss(name);
	for (int i = 0; i < size.size(); i ++) ss << size[i] << endl;
	ofstream ss2(name2);
	for (int i = 0; i < u_pop.size(); i ++) ss2 << u_pop[i] << endl;

	ss.close();
	ss2.close();
}


void save(const char* dirname, int t)
{
	char name[256],name2[256],name3[256];
	sprintf(name, "%sp%d.dat", dirname, t-1);
	ofstream file(name);

	myvec center = vec2(L/2,L/2);/* ���ĵ� */
	urban_type urban;
	for (int x = 0; x < L; x++)
	{
		for (int y = 0; y < L; y ++)
		{
			myvec p = vec2(x,y);
			file << pop[p] << '\t';
			if (pop[p])
			{
				urban[p] = t_neighbour();  
				urban[p].pop_unit = pop[p];
			}
		}
		file << endl;
	}
	file.close();

	sprintf(name, "%sz%d.dat", dirname, t-1);
	sprintf(name2, "%szp%d.dat", dirname, t-1);

	urban_size(name, name2, urban);
	int flag_temp=urban[(L-1)/2*L+(L-1)/2].flag;
	double dist_temp=0;

	t_sta sta;
	sta.init((L+1)/2,(L+1)/2);
	int count = 0;
	for (int x = 0; x < L; x++)
	{
		for (int y = 0; y < L; y ++)
		{ 
			int pindex = vec2(x, y);
			if(urban[pindex].flag==flag_temp)
			{
				sta.add(x,y, pop[pindex]);
				count ++;
			}
		}
	}

	sprintf(name, "%sr%d.dat", dirname, t-1);
	sta.save(name);

	sprintf(name, "%sap%d.dat", dirname, t-1);
	ofstream file_temp(name);
	myvec center_temp = vec2(L/2,L/2);
	urban_type urban_temp;

	for (int x = 0; x < L; x++)
	{
		for (int y = 0; y < L; y ++)
		{
			myvec p = vec2(x,y);
			file_temp << avg_pop[p] << '\t';
			if (avg_pop[p])
			{
				urban_temp[p] = t_neighbour();  
				urban_temp[p].pop_unit = avg_pop[p];
			}
		}
		file_temp << endl;
	}
	file_temp.close();

	sprintf(name, "%saz%d.dat", dirname, t-1);
	sprintf(name2, "%sazp%d.dat", dirname, t-1);

	urban_size(name, name2, urban_temp);
	flag_temp=urban_temp[(L-1)/2*L+(L-1)/2].flag;
	dist_temp=0;

	t_sta sta_temp;
	sta_temp.init((L+1)/2,(L+1)/2);
	count = 0;
	for (int x = 0; x < L; x++)
	{
		for (int y = 0; y < L; y ++)
		{ 
			int pindex = vec2(x, y);
			if(urban_temp[pindex].flag==flag_temp)
			{
				sta_temp.add(x,y, avg_pop[pindex]);
				count ++;
			}
		}
	}

	sprintf(name, "%sar%d.dat", dirname, t-1);
	sta_temp.save(name);

	for (int i = 0; i < A; i++)
	{
		avg_pop[i] = 0; 
	}

}


void update(int t)
{
	int cur = t%window_history;
	pa_time.clear();
	ex_time.clear();
	pa_num.clear();
	ex_num.clear();
	
	for (map<myvec, AliasTable>::iterator j = tables.begin(); j!= tables.end();  ++ j)
	{
		j->second.deletTable();
	}
	tables.clear();

	for (int i = 0; i < A; i++)
	{
		pop[i] += pop_temp[i].combine(plus<int>()); 
		pop_temp[i].clear();  
		snap_shots[cur][i] = pop[i];
		pre_pop[i] = 0;
		for(int j=0; j<window_history; j++)
			pre_pop[i] += snap_shots[j][i];
		
		if(pre_pop[i]>pop_threshold)
		{
			tables.insert(map<myvec, AliasTable>::value_type(i, AliasTable()));
		}	
	}

	
	for (map<myvec, AliasTable>::iterator j = tables.begin(); j!= tables.end();  ++ j)
	{
		myvec pos = j->first;
		double *prob;
		int *tag;
		tag = (int *) malloc(candidate_size*sizeof(int));
		prob = (double *) malloc(candidate_size*sizeof(double));
		int tempx = pos%L;
		int tempy = pos/L;
		int x = 0;
		int y = 0;
		int sum = 0;

		for(int i = 0;i<cutoffRegion.cur;i++)
		{
			x = tempx + cutoffRegion.x[i];
			y = tempy + cutoffRegion.y[i];
			if(x>=L)
				x -= L;
			else if(x<0)
				x += L;

			if(y>=L)
				y -= L;
			else if(y<0)
				y += L;

			tag[i] = vec2(x, y);
	
			prob[i] = (pre_pop[vec2(x, y)]+rho0)*pow(cutoffRegion.dist[i], alpha);
			sum += prob[i];

		}
		
		if(sum==0)
		{
			for(int i = 0;i<cutoffRegion.cur;i++)
			{
				x = tempx + cutoffRegion.x[i];
				y = tempy + cutoffRegion.y[i];

				if(x>=L)
					x -= L;
				else if(x<0)
					x += L;

				if(y>=L)
					y -= L;
				else if(y<0)
					y += L;

				prob[i] = pow(cutoffRegion.dist[i], alpha);
				tag[i] = vec2(x, y);
			}				
		}
		
		j->second.InitAliasTable(cutoffRegion.cur, prob, tag);
		free(prob);
		free(tag);
	}

}

int run(const char* dirname)
{

	t_user *user;
	
	for (int i = 0; i < window_history; i ++) 
	{
		for (int j=0; j<A; j ++)
		{
			snap_shots[i][j] = 0;
		}
	}

	user = (t_user *)malloc(N*sizeof(t_user));
	if (user==NULL)
		cout<<"Fail to allocate memory for users!"<<endl;

	unsigned long *seed_drand;
	seed_drand = (unsigned long *)malloc(N*sizeof(unsigned long));
	double pt, et;
	int pn, en;
	/////////////////////////////////
	clock_t start,temp1,temp2,temp3,finish;
	double totaltime;
	/////////////////////////////////

	/////////////////////////////////***************
	fill_n(pop, A, 0);
	fill_n(pre_pop, A, 0);
	fill_n(avg_pop, A, 0);
	initCutoffRegion();

	srand48((long)time(0));///////////////////////////////////
	for (int i = 0; i < N; i++) 
	{
		seed_drand[i] = lrand48();  
	}

	for (int i = 0; i < N; i ++) 
	{
		user[i].srand48((long)seed_drand[i]);
		user[i].init(vec2(L/2,L/2)); 
	}

	/////////////////////////////////***************

	start=clock();
	for (int t = 1; t <= 20001; t ++)  
	{
		temp1=clock();
		update(t); 
		temp2=clock();

		parallel_for ( 0, N, [&](int i)  
		{
			user[i].run_one_step(t);
		});
		temp3=clock();

		pt = pa_time.combine(plus<double>());
		et = ex_time.combine(plus<double>());
		pn = pa_num.combine(plus<int>());
		en = ex_num.combine(plus<int>());

		cout << t <<'\t'<<(double)(temp2-temp1)/CLOCKS_PER_SEC<<'\t'<<(double)(temp3-temp2)/CLOCKS_PER_SEC<< "\t( "<<et<<", "<<en<<")\t"<<"( "<<pt<<", "<<pn<<")"<<endl;

		if((t %500)>=(500-window_output))
		{
			for (int i = 0; i < A; i++)
			{
				avg_pop[i] += pop[i];
			}
		}
		if (t %500 == 1&& t >1)
		{
			save(dirname, t);
		}	

	}

	pa_time.clear();
	ex_time.clear();
	pa_num.clear();
	ex_num.clear();
	for (int i = 0; i < A; i++)
	{
		pop_temp[i].clear();
	}

	delete []user;
	delete []seed_drand;
}

int main(int argc, char* argv[])
{	
	snap_shots = new int*[window_history];
	int default_N = 30000;
	double default_rho0 = 0.000000000001;
	double default_alpha = -2.55;

	for (int i = 0; i < window_history; i ++) 
	{
		snap_shots[i] = new int[A];
	}

	memory_flag =1;

	N = default_N;
	alpha = default_alpha;
	rho0 = default_rho0;

	for (int i = 0; i < 5; i ++) 
	{
		char name[256];
		sprintf(name, "%d_%d_%.2f_%f_%d_%d\\", i, N,abs(alpha), rho0, memory_flag, window_output);		
		mkdir(name);	
		cout<<name<<"\n";
		run(name);
	}

	delete []snap_shots;
	return 0;
}
