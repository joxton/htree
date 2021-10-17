#include <R.h>
#include <math.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "uthash.h"
#include "utarray.h"



#define EPSILON 0.000001
#define MAXCAT 50
#define BIG_NUMBER 999999

// structures for fast historical splitting
struct nid{
  int id;
  double y;
  double time;
  UT_hash_handle hh;
};


// --- structures for historical variables (summaries of)
struct histvar{

  int xindx; // index in historical variable, ie  daop.x[xindx][] is the variable
  double *delta;
  double *delta0;
  int ndelta;
  
  double *quantiles;  // vector of probabilities [0,1]
  int nquantiles;

  
  int mean; // if 1 then mean summary
  int ecdf; // if 1 then split based on ecdf


  double **s; // the matrix of summaries 
  int **nobs;
  int **start_vec;
  int *max_nobs;
  double **deltav; // lag time record
};





// structures for keeping track of node affiliation
struct rnumber{
  int row_number;
  UT_hash_handle hh;
}; 

struct node{
  int node_number;
  int nrows; // number of rows stored in hash
  UT_array *rows_array;
  struct rnumber *rows;
  UT_hash_handle hh;
};

struct nodelist{ // used to store multiple node structures
  int list_number;
  struct node *nodes;
  UT_hash_handle hh;

};

struct node *nodes=NULL;
struct node *oob=NULL; // for out-of-bag 
struct node *id_rows=NULL; // row numbers associated with each id


struct nodelist *tnode_row=NULL; // Stores rows corresponding each terminal node for each tree, used for quantile estimation



// ---------------------------------------------------------



struct split_result{
  // structure for packing results of a split 
  int success; // 0 if unsuccessful, else 1
  int varindx;
  double nL;
  double nR;		
  double predL;
  double predR;
  double point;
  double error;
  double error_reduction;
  double tau;
  double delta;
  double delta0; // used in window summary (ie prev obstime within [delta0,delta] current obstime)
  //
  double sumL;
  double sumR;
  double sumsqL;
  double sumsqR;
  int missing; // used with "missing together" approach    
               // also used for indicating whether split is done
               // on number of observations or fraction of obs above/below
               // some threshold (missing=1, =0 respectively)

  int scolumn;
  int scolumn2;
  int hsindx;
  
  // CLASSIFICATION
  double sumL_cat[MAXCAT];
  double sumR_cat[MAXCAT];
  
  
};
void predict_basic(int *ntrees,double *pred);
double tnode_basic(int tree_indx,int row_number);
double get_summary_basic(int i,int j,int m,double q);
int find_tnode(int tree_indx,int row_number);
int  sample_orderstat(int *v,int n);
void quantile_window(struct histvar *h,int dindx);
void print_concurrent();
void permute(int *p,int n);
void permute_d(double *p,int n);
void delta_window(int nodeid,double *delta,double *delta0);
void assign_method(int *method); 
int majority_vote(double *x);



void generate_summary_mean();
void mean_window(struct histvar *h,int dindx,int oc);
void basic_split_presort(double *x,int *sindx,double *y,int n,int *train_indicator,int *xi,int minnodesize,struct split_result *s);
void predict_basic_all(int *ntrees,double *pred);
void predict_basic_all_gini(int *ntrees,double *pred);
void predict_basic_gini(int *ntrees,double *pred);
void update_terminal_nodes_wrapper(int treeindx);
void update_terminal_nodes(int treeindx);
double node_prediction_wrapper(int nodeid);
double node_prediction_logistic(int nodeid);
void tnode_basic_gini(int tree_indx,int row_number,double *pred);



struct data_options{

  double **split_matrix;

  struct histvar *hs;
  double *delta;
  double *delta0;
  int ndelta;
  int dtry;
  int qtry;
  int type;
  int *delta_sample; 
  int s_allocated;
  int nquantiles;
  double *quantiles;
  double *quantile_sample;
  int *order_sample;
  int *oob_tree;     // oob for each tree
  double *pred_tree; // prediction from each tree
  int *tree_offset;
  int family;
  
  //
  int nsamp_old;

  // delta
  int set_delta;
  double *delta_unique;
  int window_delta; // regulates how (delta0,delta)-candidates are generated
  int window_summary; // 1 if window summary, else 0 
  int n_max;  // grid length for tilde splitting
  // total loop count 
  int ncc;

  int ncat; // number of categories if classification
  
  double **x;
  double *response;
  int n;
  int p;
  int yindx;
  double *time;
  int *id;

  double *lower;
  double *upper;

  double lower_time;
  double upper_time;

  // scratch space
  double *daux1,*daux2,*daux3,*daux4,*daux5,*daux6;
  int *iaux1,*iaux2,*iaux3,*iaux4,*iaux5,*iaux6,*iaux7;  

 

  
  // predictions
  double *predictions;
  double *predictions_response;
  int *oob;
  double **predictions_gini;



  double mean_y;
  // split info and parameters 
  int min_nodesize;
  int nsamp;
  int nsamp_window;  // for window method .. 
  int nsamp_var;
  double lambda;
  double best_sumsq;
  double best_frac;
  double best_cut;
  double best_delta;
  int best_n_right;
  int best_n_left;
  int vtype;
  int valid_split;
  int subsample;
  double sample_fraction;
  int random_split;
  int method;
  int method_family;
  
  int *node;
  int nboost;
  int nsplit;
  int rf;  // if rf=1, then random forest, else boosting
  int mtry;
  int time_split; // determines whether splitting is done on time ... if =1, then yes, otherwise no 
  double delta_multiple; // multiply delta   
  // oob matrix
  //int **oob;
  int max_obs_history;
  
  // Tree matrix
  double **tree_matrix;
  int row_counter;
  int tree_counter;
  int *tree_start;
  int *tree_end;

  int *train;

  double **varimp; // varimp[k][j] is variable importance of k-th variable on j-th tree. 

  // information regarding historical and concurrent variable status
  int *splitvar_history;
  int *splitvar_concurrent;	
  int nvar_history;
  int nvar_concurrent;
  double **n_row_number;
  double **n_row_accumulator;
  int **n_changed;       // keeps track of which n_i(tau) change (for given tau)
  int *counter_changed; // keeps track of number of n_i(tau) changed (given tau)

  struct split_result **sr; 


  // mean/sum historical summaries
  double **xmat;
  double **nxmat;
  int meansummary;
  int *column_double;
  double *y_meansum;
  double LARGE_NUMBER; // Used for "Missing Together" splitting 
  
};

struct data_options daop;

void sampleWOR(int n,int size,int *res)
{
  int i,k,j,nElements;
  int *sampVec=(int*) malloc(n*sizeof(int));


  double nA,h;

  for(i=1;i<=n;i++)
    sampVec[i-1]=i-1;

  nElements=n;
  nA=(double) n;
 
  for(i=0;i<size;i++)
    {

      h=nA*unif_rand(); 
#ifdef DEBUG
      Rprintf("sampleWOR-function: unif=%lf \n",h);
#endif /*DEBUG*/

      k=(int) h;
      res[i]=sampVec[k];

      for(j=(k+1);j<nElements;j++)
	{
	  sampVec[j-1]=sampVec[j];
	}
      nElements=nElements-1;
      nA=nA-1;
    }

  free(sampVec);

}






double test_error()
{
  double error,h;
  int i,ntest;
  ntest=0;
  error=0;

  if(daop.rf!=1)
    {
  for(i=0;i<daop.n;i++)
    {
      if(daop.train[i]==0)
	{
	  error+=(daop.response[i]*daop.response[i]);
	  ntest++;
	}
    }
   }else{
    
  for(i=0;i<daop.n;i++)
    {
      if(daop.oob[i]>0)
	{
	  h=daop.response[i]-daop.predictions[i]/((double) daop.oob[i]);
	  error+=(h*h);
	  ntest++;
	}
    }
  }
  if(ntest>0)
    error=error/((double) ntest);

  return(error);
}




double train_error()
{
  double error,h;
  double pr,y;
  int i,ntest;
  ntest=0;
  error=0;

  if(daop.rf!=1)
    {
  for(i=0;i<daop.n;i++)
    {
      if(daop.train[i]==1)
	{
	  if(daop.family==1)
	    error+=(daop.response[i]*daop.response[i]);
	  if(daop.family==2)
	    {
	      pr=1/(1+exp(-daop.predictions[i]));
	      y=daop.x[daop.yindx][i];
	      error+=(-y*log(pr)-(1-y)*log(1-pr));
	    }
	  ntest++;
	}
    }
   }else{
    
  for(i=0;i<daop.n;i++)
    {
      if(daop.oob[i]>0)
	{
	  h=daop.response[i]-daop.predictions[i]/((double) daop.oob[i]);
	  error+=(h*h);
	  ntest++;
	}
    }
  }
  if(ntest>0)
    error=error/((double) ntest);

  return(error);
}




double test_error_gini()
{
  double error,h;
  int i,ntest,yint,hp;
  ntest=0;
  error=0;


    
  for(i=0;i<daop.n;i++)
    {
      if(daop.oob[i]>0)
	{
	  hp=majority_vote(daop.predictions_gini[i]);
	  h=1;
	  yint=(int) daop.response[i];
	  if(yint==hp)
	    h=0;
	  // h=daop.response[i]-daop.predictions[i]/((double) daop.oob[i]);
	  error+=h; // (h*h);
	  ntest++;
	}
    }
  
  if(ntest>0)
    error=error/((double) ntest);

  return(error);
}


int majority_vote(double *x)
{
  int class,k;
  double mv;
  
  class=0;
  mv=0;
  for(k=0;k<daop.ncat;k++)
    {
      if(mv<x[k])
	{
	  mv=x[k];
	  class=k;
	}

    }
  return(class);

}




int find_tnode(int tree_indx,int row_number)
{
  // For subject in row_number, find terminal node for subject, return terminal node number

  int k,tnode,go,j,i,split_var;
  double ff;
  int counter;
  int pred;

  tnode=-99;
  k=daop.tree_start[tree_indx];
  go=1;
  counter=0;
   while(go==1)
    {

      if(daop.tree_matrix[k][4]<(-1))
	{
	  // terminal node
	  pred=((int) daop.tree_matrix[k][2]);
	  break;
	}else{
	split_var=((int) daop.tree_matrix[k][4]);

	if(daop.tree_matrix[k][5]<(-1)) // not a summary variable
	  {
	    if(daop.x[split_var][row_number]<daop.tree_matrix[k][6])
	      {
		// Go Left
		k=((int) daop.tree_matrix[k][8]);
	      }else{
	      // Go Right
	      k=((int) daop.tree_matrix[k][9]);
	    }

	  }else{
	 
	  if(ff<daop.tree_matrix[k][5])
	    {
		// Go Left
		k=((int) daop.tree_matrix[k][8]);
	      }else{
	      // Go Right
	      k=((int) daop.tree_matrix[k][9]);
	    }
	}

      }
    }

   return(pred);
}



void subsample()
{
  int j;
  //GetRNGstate();

 for(j=0;j<daop.n;j++)
   {
     if(unif_rand()<=daop.sample_fraction)
       {
	 daop.train[j]=1;
       }else{
       daop.train[j]=0;
     }
   }

 // PutRNGstate();
}











void add_bag(int treeindx,int r_number)
{
  struct node *s;
  struct rnumber *r;

 
  HASH_FIND_INT(oob,&treeindx,s);
  if(s==NULL)
    {
      s=(struct node*) malloc(sizeof(struct node));
      s->node_number=treeindx;
      s->rows=NULL;
      s->nrows=0;
      HASH_ADD_INT(oob,node_number,s);
      //utarray_new(s->rows_array,&ut_int_icd);
    }

  HASH_FIND_INT(s->rows,&r_number,r);
  if(r==NULL)
    {
      r=(struct rnumber*) malloc(sizeof(struct rnumber));
      HASH_ADD_INT(s->rows,row_number,r);
      r->row_number=r_number;
      s->nrows++;
      //utarray_push_back(s->rows_array,&r_number); 
    }

 
}





void add_id_row(int id_indx,int r_number)
{
  struct node *s;
  struct rnumber *r;

 
  HASH_FIND_INT(id_rows,&id_indx,s);
  if(s==NULL)
    {
      s=(struct node*) malloc(sizeof(struct node));
      s->node_number=id_indx;
      s->rows=NULL;
      s->nrows=0;
      HASH_ADD_INT(id_rows,node_number,s);
      //utarray_new(s->rows_array,&ut_int_icd);
    }

  HASH_FIND_INT(s->rows,&r_number,r);
  if(r==NULL)
    {
      r=(struct rnumber*) malloc(sizeof(struct rnumber));
      HASH_ADD_INT(s->rows,row_number,r);
      r->row_number=r_number;
      s->nrows++;
      //utarray_push_back(s->rows_array,&r_number); 
    }

 
}



void get_oob(int *oob_vec)
{

  //
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  int k;

  for(k=0;k<(daop.nboost*daop.n);k++)
    oob_vec[k]=0;

  k=0;
  HASH_ITER(hh,oob,s,stemp){
    HASH_ITER(hh,s->rows,r,rtemp){
      oob_vec[k*daop.n+r->row_number]=1;
    }
	k++;
  }

}




void subsample_oob(int treeindx)
{
  int j;
  // Checked
  //GetRNGstate();

 for(j=0;j<daop.n;j++)
   {
     if(unif_rand()<=daop.sample_fraction)
       {	
	 daop.train[j]=1;
	 
       }else{
       daop.train[j]=0;
        add_bag(treeindx,j);
     }
   }


 //PutRNGstate();
}


void subsample_id(int treeindx)
{
 
  int j;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  
  //GetRNGstate();

  // for(j=0;j<daop.n;j++)
  HASH_ITER(hh,id_rows,s,stemp)
    {
      if(unif_rand()<=daop.sample_fraction)
	{
	  HASH_ITER(hh,s->rows,r,rtemp){	 
	    daop.train[r->row_number]=1;	  
	  }
	}else{
	HASH_ITER(hh,s->rows,r,rtemp){	 
	  daop.train[r->row_number]=0;
	  add_bag(treeindx,r->row_number);
	}
      
       
      }
    }


  //PutRNGstate();


}



void delete_oob()
{

  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  HASH_ITER(hh,oob,s,stemp){
    HASH_ITER(hh,s->rows,r,rtemp){
      HASH_DEL(s->rows,r);
      free(r);
    }
    HASH_DEL(oob,s);
    free(s);
  }

}


void delete_id_rows()
{

  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  HASH_ITER(hh,id_rows,s,stemp){
    HASH_ITER(hh,s->rows,r,rtemp){
      HASH_DEL(s->rows,r);
      free(r);
    }
    HASH_DEL(id_rows,s);
    free(s);
  }

}



void set_oob(int *oob_vec)
{
  int i,b;
  delete_oob();

  daop.subsample=0; // turn it off 
  for(b=0;b<daop.nboost;b++)
    {
      for(i=0;i<daop.n;i++)
	{
	  if(oob_vec[b*daop.n+i]==1)
	    add_bag(b,i);
	}

    }
  

}


//void add_id_row(int id_indx,int r_number)


void set_id_rows(int *id_vec)
{
  int i;
  delete_id_rows();

      for(i=0;i<daop.n;i++)
	{
	  add_id_row(id_vec[i],i);
	}

}



void add_row(int n_number,int r_number)
{
  struct node *s;
  struct rnumber *r;
  HASH_FIND_INT(nodes,&n_number,s);
  if(s==NULL)
    {
      s=(struct node*) malloc(sizeof(struct node));
      s->node_number=n_number;
      HASH_ADD_INT(nodes,node_number,s);
      s->rows=NULL;
      s->nrows=0;
      utarray_new(s->rows_array,&ut_int_icd);
    }

  HASH_FIND_INT(s->rows,&r_number,r);
  if(r==NULL)
    {
      r=(struct rnumber*) malloc(sizeof(struct rnumber));
      HASH_ADD_INT(s->rows,row_number,r);
      r->row_number=r_number;
      s->nrows++;
      utarray_push_back(s->rows_array,&r_number); 
    }

}

void delete_node(int n)
{
  struct node *s;
  struct rnumber *r,*rtemp;

  HASH_FIND_INT(nodes,&n,s);
  if(s==NULL){}else{
    HASH_ITER(hh,s->rows,r,rtemp)
      {
	HASH_DEL(s->rows,r);
	free(r);
      }
    utarray_free(s->rows_array);
    
    HASH_DEL(nodes,s);
    free(s);
  }
}

void delete_all_nodes()
{
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  HASH_ITER(hh,nodes,s,stemp){
      HASH_ITER(hh,s->rows,r,rtemp)
	{
	  HASH_DEL(s->rows,r);
	  free(r);
	}
      utarray_free(s->rows_array);
      HASH_DEL(nodes,s);
      free(s);
    }

}





int sample_indx(int n)
{
  // sample an index from 0:(n-1)
  int i;
  double nd,u,uu;
  u=unif_rand();
  nd=((double) n);
  uu=nd*u;
  i=((int) uu);

  if(i==n)
    i=n-1;
  
  return(i);
}

void sample_indx_wrapper(int *n,int *res)
{

 GetRNGstate();

 res[0]=sample_indx(n[0]);
 
 PutRNGstate();
  
}


void permute(int *p,int n)
{

  // Take p and permute it
  int k,i,pa;
  
  for(k=0;k<n;k++)
    {
      i=sample_indx((n-k));
      pa=p[(n-k-1)];
      p[(n-k-1)]=p[i];
      p[i]=pa;
    }

}


void permute_d(double *p,int n)
{

  // Take p and permute it (p-double vector)
  int k,i;
  double pa;
  
  for(k=0;k<n;k++)
    {
      i=sample_indx((n-k));
      pa=p[(n-k-1)];
      p[(n-k-1)]=p[i];
      p[i]=pa;
    }

}



void permute_wrapper(int *n,int *res)
{
  int *p;
  int k;
  //p=(int*) malloc(sizeof(int)*n[0]);
  for(k=0;k<n[0];k++)
    res[k]=k;
  
  GetRNGstate();
  permute(res,n[0]);
  PutRNGstate();

}


double node_prediction_uthash(int nodeid)
{
  // return node mean
  int i;
  double a;
  double n;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  
  a=0;
  n=0;

  HASH_FIND_INT(nodes,&nodeid,s);
  if(s==NULL){}else{
   
    HASH_ITER(hh,s->rows,r,rtemp)
      {
	i=r->row_number;
	if(daop.train[i]==1)
	  {
	    a+=daop.response[i];
	    n++;
	  }
       
      }
  }
   
  if(n<.5)
    n=1;
  a=a/n;

  return(a);
}




void node_prediction_gini(int nodeid,double *a)
{
  // return node mean
  int i,yint;
  //double a;
  double n;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  for(i=0;i<daop.ncat;i++)
    a[i]=0;
  
  n=0;

  HASH_FIND_INT(nodes,&nodeid,s);
  if(s==NULL){}else{
    //for(i=0;i<daop.n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {
      i=r->row_number;
	  if(daop.train[i]==1)
	    {
	      yint=(int) daop.response[i];
	      a[yint]++; //=daop.response[i];
	      n++;
	    }
       
    }
  }


}







void set_train_oob(int treeindx)
{

  // Set daop.train
  int i,b;
  struct node *s;
  struct rnumber *r,*rtemp;
  
  for(i=0;i<daop.n;i++)
    daop.train[i]=1;


  HASH_FIND_INT(oob,&treeindx,s);
  if(s==NULL){Rprintf(" Did not find oob-list element ...\n");}else{

    HASH_ITER(hh,s->rows,r,rtemp)
      {
	daop.train[r->row_number]=0;
      }

  }

  
}



void initialize_split_result(struct split_result *spr)
{
  int k;
  spr->success=0;
  spr->nL=0;
  spr->nR=0;
  spr->sumL=0;
  spr->sumR=0;
  spr->sumsqL=0;
  spr->sumsqR=0;
  spr->predL=0;
  spr->predR=0;
  spr->point=0;
  spr->error=0;
  spr->delta0=-99;
  spr->missing=0;
  spr->delta=-99;
  spr->tau=-99;
  spr->scolumn=-99;
  spr->scolumn2=-99;
  spr->hsindx=-99;
  
  if(daop.ncat>0)
    {
      // ncat>0 => classification
      for(k=0;k<daop.ncat;k++)
	{
	  spr->sumL_cat[k]=0;
	  spr->sumR_cat[k]=0;
	}

    }
}

void error_split_result(struct split_result *spr)
{

  // Compute error of split result (also predL and predR are formed) 
  spr->predL=spr->sumL/spr->nL;
  spr->predR=spr->sumR/spr->nR;
  spr->error=(spr->nL/(spr->nL+spr->nR))*(spr->sumsqL/spr->nL-spr->predL*spr->predL)+(spr->nR/(spr->nL+spr->nR))*(spr->sumsqR/spr->nR-spr->predR*spr->predR);



}



void error_split_gini(struct split_result *spr)
{

  // Compute error of split result (also predL and predR are formed)
  double eL,eR,p;
  int k;
  eL=0;
  eR=0;
  for(k=0;k<daop.ncat;k++)
    {
      p=spr->sumL_cat[k]/spr->nL;
      eL+=(p*(1-p));
      p=spr->sumR_cat[k]/spr->nR;
      eR+=(p*(1-p));
      

    }
  spr->predL=eL;
  spr->predR=eR;
  spr->error=(spr->nL/(spr->nL+spr->nR))*eL+(spr->nR/(spr->nL+spr->nR))*eR;


}


void error_reduction_split_result(struct split_result *spr)
{
  // Compute error reduction, this is error prior to node splitting minus error post split 
  double error_prior;

  error_prior=(spr->sumL+spr->sumR)/(spr->nL+spr->nR);
  error_prior=error_prior*error_prior;
  error_prior=(spr->sumsqL+spr->sumsqR)/(spr->nL+spr->nR)-error_prior;
  
  spr->error_reduction=error_prior-spr->error;
      
}

void error_reduction_split_gini(struct split_result *spr)
{
  // Compute error reduction, this is error prior to node splitting minus error post split 
  double error_prior;
  int k;
  double p;

  error_prior=0;
  for(k=0;k<daop.ncat;k++)
    {
      p=(spr->sumL_cat[k]+spr->sumR_cat[k])/(spr->nL+spr->nR);
      error_prior+=(p*(1-p));

    }

  spr->error_reduction=error_prior-spr->error;
      
}





void basic_split(double *x,double *y,int n,int *xi,int minnodesize,struct split_result *s)
{
  // basic splitting. given predictor x, response y, find best split point, return in structure split_result

  // Splitting is defined as:
  // x<best_point -> go left
  // x>=best_point -> go right 
  
  double sumL,sumR,nL,nR,sumsqL,sumsqR,predL,predR,split_error;
  int k;
  double min_error,best_nL,best_nR,best_predL,best_predR,best_point,best_ssqL,best_ssqR,best_sumL,best_sumR;
  double xc,yaux,mnsd;
  int init;
  struct split_result csplit,bsplit;

 initialize_split_result(&csplit);
  
  mnsd=(double) minnodesize;
  
  for(k=0;k<n;k++){
    xi[k]=k;
    csplit.nL++;
    csplit.sumL+=y[k];
    csplit.sumsqL+=(y[k]*y[k]);
  }
  
  revsort(x,xi,n);  // sorts descending order

  xc=x[0];
  csplit.point=xc;
  init=0;
  for(k=0;k<n;k++)
    {

      if(fabs(x[k]-xc)>EPSILON)
	{

	  if(csplit.nL>=mnsd&&csplit.nR>=mnsd)
	    {
	      // valid split
	      error_split_result(&csplit); // compute error at split point 	
	      if(init==0||(csplit.error<bsplit.error))
		{
		  bsplit=csplit;
		  init=1;
		}
	    }
	  
	}
      csplit.nR++;
      csplit.nL--;
      yaux=y[xi[k]];
      csplit.sumL-=yaux;
      csplit.sumR+=yaux;
      csplit.sumsqR+=(yaux*yaux);
      csplit.sumsqL-=(yaux*yaux);
      xc=x[k];
      csplit.point=xc;
      
      if(csplit.nL<mnsd)
	break;
    }


  if(init!=0)
    {
      s[0]=bsplit;
      s->success=1;
      
    }else{
    s->success=0;
    }
  
}





void basic_split_gini(double *x,double *y,int n,int *xi,int minnodesize,struct split_result *s)
{
  // basic splitting. given predictor x, response y, find best split point, return in structure split_result

  // Splitting is defined as:
  // x<best_point -> go left
  // x>=best_point -> go right 
  
  double sumL,sumR,nL,nR,sumsqL,sumsqR,predL,predR,split_error;
  int k;
  double min_error,best_nL,best_nR,best_predL,best_predR,best_point,best_ssqL,best_ssqR,best_sumL,best_sumR;
  double xc,yaux,mnsd;
  int init;
  int yint;
  struct split_result csplit,bsplit;

  initialize_split_result(&csplit);  // should work for classification
  
  mnsd=(double) minnodesize;
  
  for(k=0;k<n;k++){
    xi[k]=k;
    csplit.nL++;
    yint=(int) y[k];
    csplit.sumL_cat[yint]++; //+=y[k];
    //csplit.sumsqL+=(y[k]*y[k]);
  }
  
  revsort(x,xi,n);  // sorts descending order

  xc=x[0];
  csplit.point=xc;
  init=0;
  for(k=0;k<n;k++)
    {

      if(fabs(x[k]-xc)>EPSILON)
	{

	  if(csplit.nL>=mnsd&&csplit.nR>=mnsd)
	    {
	      // valid split
	      error_split_gini(&csplit); // compute error at split point
	      if(init==0||(csplit.error<bsplit.error))
		{
		  bsplit=csplit;
		  init=1;
		}
	    }
	  
	}
      csplit.nR++;
      csplit.nL--;
      yint=(int) y[xi[k]];
      csplit.sumL_cat[yint]--; //=yaux;
      csplit.sumR_cat[yint]++; //=yaux;
      xc=x[k];
      csplit.point=xc;
      
      if(csplit.nL<=mnsd)
	break;
    }


  if(init!=0)
    {
      s[0]=bsplit;
      s->success=1;
      
    }else{
    s->success=0;
    }
  
}









void split_standard(int nodeid,struct split_result *best_split)
{
  // split node in standard manner

  struct split_result sr;
  double *x,*y;
  int *xi;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j,mtry_use;
  int *predictors;
  predictors=(int*) malloc(sizeof(int)*daop.p*2);

  
  x=daop.daux5;
  y=daop.daux4;
  xi=daop.iaux5;
  init=0;
  best_split->success=0; // until proven otherwise....


  counter=0;
  for(k=0;k<daop.nvar_concurrent;k++)
    {
      if(daop.splitvar_concurrent[k]!=daop.yindx)
	{

	  predictors[counter]=daop.splitvar_concurrent[k];
	  counter++;
	}
    }
  permute(predictors,(counter));
  
  mtry_use=daop.mtry;
  if(mtry_use>daop.nvar_concurrent)
    mtry_use=daop.nvar_concurrent;
  
  for(j=0;j<mtry_use;j++)
    {
	k=predictors[j];
      

      if(k!=daop.yindx)
	{
	  // get data
	  HASH_FIND_INT(nodes,&nodeid,s);

	  counter=0;
	  if(s==NULL){}else{ 

	    HASH_ITER(hh,s->rows,r,rtemp)
	      {
		i=r->row_number;
		if(daop.train[i]==1){
		  //y[counter]=daop.x[daop.yindx][i];
		  y[counter]=daop.response[i];
		  x[counter]=daop.x[k][i];
		  counter++;
		}
	      }
	  }
	  if(counter>0)
	    {
	      basic_split(x,y,counter,xi,daop.min_nodesize,&sr);
	      if(sr.success==1)
		{
		  if(init==0||((best_split->error)>sr.error))
		    {
		      best_split[0]=sr;
		      best_split->varindx=k;
		      init=1;
		    }
		  
		}
	      

	    }
	}
      
    }

  free(predictors);
}









void split_standard_gini(int nodeid,struct split_result *best_split)
{
  // split node in standard manner

  struct split_result sr;
  double *x,*y;
  int *xi;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j,mtry_use;
  // int predictors[1000];
  int yint;
  int *predictors;

  predictors=(int*) malloc(sizeof(int)*daop.n);
  
  x=daop.daux5;
  y=daop.daux4;
  xi=daop.iaux5;

  init=0;
  best_split->success=0; // until proven otherwise....

  counter=0;
  for(k=0;k<daop.nvar_concurrent;k++)
    {

      if(daop.splitvar_concurrent[k]!=daop.yindx){
	  predictors[counter]=daop.splitvar_concurrent[k];
	  counter++;
      }
    }

  permute(predictors,counter);
  
  mtry_use=daop.mtry;
  if(mtry_use>daop.nvar_concurrent)
    mtry_use=daop.nvar_concurrent;
  
  for(j=0;j<mtry_use;j++)
    {
	k=predictors[j];
      
      
      if(k!=daop.yindx)
	{
	  // get data
	  HASH_FIND_INT(nodes,&nodeid,s);

	  counter=0;
	  if(s==NULL){}else{ 

	    HASH_ITER(hh,s->rows,r,rtemp)
	      {
		i=r->row_number;
		if(daop.train[i]==1){
		  y[counter]=daop.response[i];
		  x[counter]=daop.x[k][i];
		  counter++;
		}
	      }
	  }
	  
	  if(counter>0)
	    {
	      basic_split_gini(x,y,counter,xi,daop.min_nodesize,&sr);
	      
	      if(sr.success==1)
		{
		  if(init==0||((best_split->error)>sr.error))
		    {
		      best_split[0]=sr;
		      best_split->varindx=k;
		      init=1;
		    }
		  
		}
	      

	    }
	}
      
    }

  free(predictors);
}








void initialize_prediction_response()
{
  int i;

  if(daop.method_family==2) // logistic 
    {

      for(i=0;i<daop.n;i++)
	{
	  daop.predictions[i]=0;
	  daop.predictions_response[i]=1/(1+exp(-daop.predictions[i]));
	  daop.response[i]=(daop.x[daop.yindx][i]-daop.predictions_response[i]);
	}

    }
  
}

double test_error_wrapper()
{

  double terr;

  terr=0;
  if(daop.method_family==1){
  if(daop.method==1||((daop.method>2)&&(daop.method<9)))
    {
      // daop.method=1: regression standard
      terr=test_error();
    }
  }
  if(daop.method_family==2)
    {
      //      terr=test_error_logistic();
    }
  // for other methods .... 

  if(daop.method==9)
    terr=test_error_gini();

  
  return(terr);
}









double test_error_logistic()
{

  
  double error,h;
  int i,ntest;
  ntest=0;
  error=0;


  for(i=0;i<daop.n;i++)
    {
      if(daop.train[i]==0)
	{
	  h=daop.x[daop.yindx][i]*log(daop.predictions_response[i]);
	  h=h+(1-daop.x[daop.yindx][i])*log((1-daop.predictions_response[i]));
	  error+=h;
	  ntest++;
	}
    }
    
  if(ntest>0)
    error=error/((double) ntest);

  return(error);
}





void update_residual_logistic(int tree_indx)
{
  int k,i,start,go,ni,ti,end;
  //double node_predictions[10000]; // no more than 100 terminal nodes
  //int tnodes[10000];// terminal nodes
  int ntnodes;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  int *tnodes;
  double *node_predictions;
  
  tnodes=(int*) malloc(sizeof(int)*daop.n);
  node_predictions=(double*) malloc(sizeof(double)*daop.n);
  
  ntnodes=0;
  start=daop.tree_start[tree_indx];
  end=daop.tree_end[tree_indx];
  for(k=start;k<=end;k++)
    {

      ti=(int) daop.tree_matrix[k][0];
      if(ti==tree_indx)
	{
	  ni=((int) daop.tree_matrix[k][2]);
	  node_predictions[ni]=daop.tree_matrix[k][3];
	  if(daop.tree_matrix[k][4]<(-1))
	    {
	      tnodes[ntnodes]=ni;
	      ntnodes++;
	    }
	  
	}else{
	go=0;
	break;
      }

    }
 
  for(k=0;k<ntnodes;k++)
    {

      ni=tnodes[k];
      HASH_FIND_INT(nodes,&ni,s);
      if(s==NULL){}else{

	HASH_ITER(hh,s->rows,r,rtemp){
    
 	  i=r->row_number;
	  daop.predictions[i]+=(node_predictions[ni]);
	  daop.predictions_response[i]=1/(1+exp(-daop.predictions[i]));
	  daop.response[i]=(daop.x[daop.yindx][i]-daop.predictions_response[i]);

	}
      }
    }

 

  free(tnodes);
  free(node_predictions);
  
}






// ----- quantile regression functions ---------------------------------



// tnode_row

void add_tnode_row(int treeindx,int nodeindx,int rowindx)
{
 


  struct nodelist *t;
  struct node *s;
  struct rnumber *r;

  HASH_FIND_INT(tnode_row,&treeindx,t);
  if(t==NULL)
    {
      t=(struct nodelist*) malloc(sizeof(struct nodelist));
      t->list_number=treeindx;
      t->nodes=NULL;
      HASH_ADD_INT(tnode_row,list_number,t);
      //utarray_new(s->rows_array,&ut_int_icd);
    }
  
 
  HASH_FIND_INT(t->nodes,&nodeindx,s);
  if(s==NULL)
    {
      s=(struct node*) malloc(sizeof(struct node));
      s->node_number=nodeindx;
      s->rows=NULL;
      s->nrows=0;
      HASH_ADD_INT(t->nodes,node_number,s);
      //utarray_new(s->rows_array,&ut_int_icd);
    }

  HASH_FIND_INT(s->rows,&rowindx,r);
  if(r==NULL)
    {
      r=(struct rnumber*) malloc(sizeof(struct rnumber));
      HASH_ADD_INT(s->rows,row_number,r);
      r->row_number=rowindx;
      
      //utarray_push_back(s->rows_array,&r_number); 
    }

 
}



void delete_tnode_row()
{
  struct nodelist *t,*ttemp;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  HASH_ITER(hh,tnode_row,t,ttemp){
  HASH_ITER(hh,t->nodes,s,stemp){
    HASH_ITER(hh,s->rows,r,rtemp){
      HASH_DEL(s->rows,r);
      free(r);
    }
   
    HASH_DEL(t->nodes,s);
    free(s);
  }
  HASH_DEL(tnode_row,t);
  free(t);
  }
}






void quantile_aux(double *prob,double *y,int *n,double *weights,double *quant)
{
  double cc,cc_new;
  int k;

  quant[0]=-99;
  cc=0;
  for(k=0;k<n[0];k++)
    {
      cc_new=cc+weights[k];
      if(cc<=prob[0]&&cc_new>=prob[0])
	{
	  quant[0]=y[k-1]+(y[k]-y[k-1])*(prob[0]-cc)/(cc_new-cc);
	  break;
	}
	    cc=cc_new;
    }
  
}



void quantile_R(double *prob,int *np,double *y,int *n,double *weights,int *m,double *quant)
{

  // find quantiles prob[1:np] for m observations, weights in weights
  int k,j,counter,i;
  double q;

  counter=0;
  for(k=0;k<m[0];k++)
    {

      for(j=0;j<np[0];j++)
	{

	  quantile_aux(&(prob[j]),y,n,&(weights[k*n[0]]),&q);
	  quant[counter]=q;
	  counter++;
	}


    }


}



void set_tnode_row()
{
  // Fill up tnode_row list, mapping rows to terminal nodes of each tree. 
  
  int i,k,tnode;
  delete_tnode_row();

  for(k=0;k<daop.nboost;k++)
    {
      for(i=0;i<daop.n;i++){
	//	tnode=find_tnode(k,i);
	add_tnode_row(k,tnode,i);
      }

    }

}

void get_rows(int *treeindx,int *nodeindx,int *res,int *n)
{
  int i,k,tnode;

  struct nodelist *t,*ttemp;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  int counter;

  HASH_FIND_INT(tnode_row,treeindx,t);
  counter=0;
  if(t==NULL){}else{

    HASH_FIND_INT(t->nodes,nodeindx,s);
    if(s==NULL){}else{
      
	      HASH_ITER(hh,s->rows,r,rtemp){
		res[counter]=r->row_number;
		counter++;
	      }
    }

  }

  n[0]=counter;
}

void get_weights(int *n,int *res)
{
  int i,k,tnode;

  struct nodelist *t,*ttemp;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;


  for(i=0;i<(n[0]*daop.n);i++)
    res[i]=0;
  
  for(i=0;i<daop.n;i++)
    {

      for(k=0;k<daop.nboost;k++)
	{

	  tnode=find_tnode(k,i);
	  HASH_FIND_INT(tnode_row,&k,t);
	  if(t==NULL){}else{
	    HASH_FIND_INT(t->nodes,&tnode,s);
	    if(s==NULL){}else{
	      HASH_ITER(hh,s->rows,r,rtemp){
		res[i*n[0]+r->row_number]+=1;
	      }
	    }
	  }
	}

    }

}







void generate_summary_mean()
{
  int k,j;
  double *vv;

  // allocate hs[k].s[]
  for(k=0;k<=daop.nvar_history;k++)
    {
     
      for(j=0;j<daop.ndelta;j++){
	 vv=(double*) malloc(sizeof(double)*(daop.n));
	daop.hs[k].s[j]=vv;
      }
    }
  daop.s_allocated=1;
  
  for(k=0;k<daop.nvar_history;k++)
    for(j=0;j<(daop.hs[k].ndelta);j++)
      mean_window(&(daop.hs[k]),j,0);

  // generate observation count features
  k=daop.nvar_history;
  for(j=0;j<(daop.hs[k].ndelta);j++)
    mean_window(&(daop.hs[k]),j,1);



 }


void regenerate_history_summary()
{
  int k,j;

  // for varimp,
  if(daop.type==1||daop.type==4)
    {
 for(k=0;k<daop.nvar_history;k++)
    for(j=0;j<(daop.hs[k].ndelta);j++)
      mean_window(&(daop.hs[k]),j,0);

  // generate observation count features
  k=daop.nvar_history;
  for(j=0;j<(daop.hs[k].ndelta);j++)
    mean_window(&(daop.hs[k]),j,1);


    }else{

  // allocate hs[k].s[]
  for(k=0;k<daop.nvar_history;k++)
    {
     
      for(j=0;j<daop.ndelta;j++){

	  
	   free(daop.hs[k].s[j]);
	   free(daop.hs[k].deltav[j]);
	   free(daop.hs[k].nobs[j]);
	   free(daop.hs[k].start_vec[j]);
 	quantile_window(&(daop.hs[k]),j);
      }
    }


  }

}


void generate_summary_quantile()
{
  int k,j;
  double *vv;

  // allocate hs[k].s[]
  for(k=0;k<daop.nvar_history;k++)
    {
     
      for(j=0;j<daop.ndelta;j++){
	quantile_window(&(daop.hs[k]),j);
      }
    }
  daop.s_allocated=1;

 
  
    // generate observation count features
  k=daop.nvar_history;
  for(j=0;j<(daop.hs[k].ndelta);j++){
    daop.hs[k].s[j]=(double*) malloc(sizeof(double)*daop.n);
    mean_window(&(daop.hs[k]),j,1);

  }


}




void mean_window(struct histvar *h,int dindx,int obs_count)
{
  
  // Find mean of all lagged values within k-th observation window, ie all obs of variable within (delta,delta0] of observation time
  // and write it to h->s[dindx][]


  // If obs_count==0 then mean is computed, else function returns number of observations 
   
  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  double nobs;
  struct data_options *d;
  double delta,delta0;
  int k,vindx;
  
  d=&daop;

  delta=daop.delta[dindx];//h->delta[dindx];
  delta0=daop.delta0[dindx];//h->delta0[dindx];

  vindx=h->xindx;
 

  for(i=0;i<d->n;i++)
    {

       
      cid=d->id[i];

      if(obs_count==0)
	h->s[dindx][i]=0;
      nobs=0;
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if(((d->time[i]-d->time[i-j])<=delta)&&((d->time[i]-d->time[i-j])>=delta0))  // used to be j-i 
	    {
	      if(obs_count==0)
		h->s[dindx][i]+=d->x[vindx][i-j];
	      nobs+=1;
	    
	    }else{
	    if((d->time[i]-d->time[i-j])>delta)
	      break;
	  }
	    
	}
      
      if(nobs>EPSILON)
	{
	  if(obs_count==0){
	    h->s[dindx][i]=h->s[dindx][i]/nobs;
          }else{
	    h->s[dindx][i]=nobs;

	  }


	}else{
	h->s[dindx][i]=0;
      }

 

    }




}





void quantile_window(struct histvar *h,int dindx)
{
  
  // Find quantiles of all lagged values within k-th observation window, ie all obs of variable within (delta,delta0] of observation time
  // and write it to h->s[dindx][]


  // If obs_count==0 then mean is computed, else function returns number of observations 
   
  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  double nobs;
  struct data_options *d;
  double delta,delta0;
  int k,vindx;
  double *pvec;
  int *nobs_vec;
  int *start_vec;
  int max_nobs;
  int *ip;
  double *daux;
  double *delta_vec;
  
  d=&daop;

  delta=daop.delta[dindx]; //h->delta[dindx];
  delta0=daop.delta0[dindx]; //h->delta0[dindx];

  vindx=h->xindx;

  // find number of obs within window
  nobs_vec=(int*) malloc(sizeof(int)*daop.n);
  start_vec=(int*) malloc(sizeof(int)*daop.n);
  //  end_vec=(int*) malloc(sizeof(int)*daop.n);

  
  max_nobs=0;
  for(i=0;i<d->n;i++)
    {
      start_vec[i]=0;
      nobs_vec[i]=0;
      cid=d->id[i];
      for(j=1;j<d->n;j++)
	{
	  

	  if((i-j)<0)
	    break;
	  if(d->id[i-j]!=cid)
	    break;
	  

	  if(((d->time[i]-d->time[i-j])<=delta)&&((d->time[i]-d->time[i-j])>=delta0))  // used to be j-i 
	    {

	      nobs_vec[i]++;
	      if(nobs_vec[i]==1)
		start_vec[i]=i-j;

	    }else{
	    if((d->time[i]-d->time[i-j])>delta)
	      break;
	  }
  

	}

      if(max_nobs<nobs_vec[i])
	max_nobs=nobs_vec[i];
    }

   pvec=(double*) malloc(sizeof(double)*daop.n*max_nobs);
  delta_vec=(double*) malloc(sizeof(double)*daop.n*max_nobs);
  
  counter=0;
  for(i=0;i<d->n;i++)
    {

       
      cid=d->id[i];
      start_vec[i]=0;
      nobs_vec[i]=0;
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if(((d->time[i]-d->time[i-j])<=delta)&&((d->time[i]-d->time[i-j])>=delta0))  
	    {
	      nobs_vec[i]++;
	      if(nobs_vec[i]==1)
		start_vec[i]=counter;
	      pvec[counter]=d->x[vindx][i-j];
	      delta_vec[counter]=(d->time[i]-d->time[i-j]);
	    

	      counter++;
	    
	    }else{
	    if((d->time[i]-d->time[i-j])>delta)
	      break;
	  }
	    
	}

    }

  // ..sort history
  daux=(double*) malloc(sizeof(double)*max_nobs);
  ip=(int*) malloc(sizeof(int)*max_nobs);
  for(i=0;i<max_nobs;i++)
    ip[i]=i;
  counter=0;
  for(i=0;i<daop.n;i++){
    if(nobs_vec[i]>1){

      if(daop.type>4){
        for(j=0;j<nobs_vec[i];j++)
      	  ip[j]=j;
      }
      revsort(&(pvec[start_vec[i]]),ip,nobs_vec[i]);
      counter+=nobs_vec[i];

     if(daop.type>4){
        for(j=0;j<nobs_vec[i];j++)
     	  daux[j]=delta_vec[start_vec[i]+j];
        for(j=0;j<nobs_vec[i];j++)
     	  delta_vec[start_vec[i]+j]=daux[ip[j]];
      }
      
    }
  }

  free(ip);
  free(daux);

   h->deltav[dindx]=delta_vec;
  h->s[dindx]=pvec;
  h->nobs[dindx]=nobs_vec;
  h->start_vec[dindx]=start_vec;
  h->max_nobs[dindx]=max_nobs; // max number of observations for this delta combination 
  
 
 

}








// -------------------------------------------------------------------
// ---------------- SUMMARIES  ---------------------------------------
// -------------------------------------------------------------------

void create_history(double *delta,double *delta0,int *ndelta,int *type,double *quantiles,int *nquantiles,int *dtry)
{

  int i,j,k,m;
  double *v;
  int max_obs;

  if(daop.nvar_history>0)
    {
       daop.dtry=dtry[0];
      daop.delta_sample=(int*) malloc(sizeof(int)*(ndelta[0]+10));
      daop.ndelta=ndelta[0];
      daop.type=type[0];
      daop.delta=delta;
      daop.delta0=delta0;

      
      for(k=0;k<=daop.nvar_history;k++)
	{


	  if(k<daop.nvar_history)
	    daop.hs[k].xindx=daop.splitvar_history[k];
	  
	  if(k==daop.nvar_history)
	    daop.hs[k].xindx=-1; // split on number of obs 
	  
	  
	  daop.hs[k].ndelta=ndelta[0];

 
	 
	    daop.hs[k].s=(double**) malloc(sizeof(double*)*ndelta[0]);
	    if(type[0]==2||type[0]==3||type[0]==5||type[0]==6){
	      daop.hs[k].nobs=(int**) malloc(sizeof(int*)*ndelta[0]);
	      daop.hs[k].start_vec=(int**) malloc(sizeof(int*)*ndelta[0]);
	      daop.hs[k].max_nobs=(int*) malloc(sizeof(int)*ndelta[0]);
	      daop.hs[k].deltav=(double**) malloc(sizeof(double*)*ndelta[0]); 
	    }
	}

      if(type[0]==1||type[0]==4)
	    {
	      generate_summary_mean();

	    }


	  if(type[0]==2||type[0]==3||type[0]==5||type[0]==6)
	    {
	      
	      daop.quantiles=quantiles;
	      daop.nquantiles=nquantiles[0];
	      
	      // quantile summary 
	      generate_summary_quantile();

	      if(type[0]==3||type[0]==6)
		{
		  max_obs=0;
		  for(k=0;k<daop.ndelta;k++)
		    {
		      for(j=0;j<daop.nvar_history;j++)
			{
			  if(daop.hs[j].max_nobs[k]>max_obs)
			    max_obs=daop.hs[j].max_nobs[k];
			}
		    }
		      
		  daop.order_sample=(int*) malloc(sizeof(int)*max_obs);
 
		}


	    }


	

    }


  
}



void sample_delta_index(int *dvec)
{
  //int dvec[1000];

  int ndelta;
  int k; 
  ndelta=daop.ndelta;

  for(k=0;k<ndelta;k++)
    dvec[k]=k;

  permute(dvec,ndelta);
  
  

 
}




void sample_quantiles(double *qvec)
{
  //int dvec[1000];

  int ndelta;
  int k; 
  ndelta=daop.ndelta;

  for(k=0;k<daop.qtry;k++)
    qvec[k]=daop.quantiles[k];

  permute_d(qvec,daop.qtry);
  
}




int  sample_quantiles2(double *qvec,int maxobs)
{
  //int dvec[1000];

  int ndelta;
  int k; 
 

  if(maxobs>=daop.nquantiles)
    {
  for(k=0;k<daop.nquantiles;k++)
    qvec[k]=daop.quantiles[k];

  permute_d(qvec,daop.qtry);
    }else{

    for(k=0;k<maxobs;k++){
      qvec[k]=((double) (k+1))/((double) maxobs); //daop.quantiles[k];
    }
     permute_d(qvec,maxobs); 

  }
  k=daop.qtry;
  if(k>maxobs)
    k=maxobs;
  
  return(k);
}





int  sample_quantiles_w(double *qvec,int m,int m2) // int maxobs)
{

  int ndelta;
  int k,maxobs;


  // this is approximate
  if(m2>(-1)){
    maxobs=daop.hs[0].max_nobs[m]-daop.hs[0].max_nobs[m2]+1;
  }else{
    maxobs=daop.hs[0].max_nobs[m];
 
  }
   
  if(maxobs>=daop.nquantiles)
    {
 
      for(k=0;k<daop.nquantiles;k++){
    qvec[k]=daop.quantiles[k];

      }

  permute_d(qvec,daop.qtry);
    }else{

    for(k=0;k<maxobs;k++){
      qvec[k]=((double) (k+1))/((double) maxobs); //daop.quantiles[k];
    }
     permute_d(qvec,maxobs); 

  }
  k=daop.qtry;
  if(k>maxobs)
    k=maxobs;
  
  return(k);
}



int  sample_orderstat(int *qvec,int maxobs)
{
  //int dvec[1000];

  int ndelta;
  int k;   
 
  for(k=0;k<maxobs;k++)
    qvec[k]=k+1;

  permute(qvec,maxobs);

  k=daop.qtry;
  if(k>maxobs)
    k=maxobs;
  
  return(k);
}





int  sample_orderstat_w(int *qvec,int m,int m2)
{
  //int dvec[1000];

  int ndelta;
  int k;
  int maxobs;

    if(m2>(-1)){
    maxobs=daop.hs[0].max_nobs[m]-daop.hs[0].max_nobs[m2]+1;
  }else{
    maxobs=daop.hs[0].max_nobs[m];
 
  }

 
  for(k=0;k<maxobs;k++)
    qvec[k]=k+1;

  permute(qvec,maxobs);

  k=daop.qtry;
  if(k>maxobs)
    k=maxobs;
  
  return(k);
}



int quantile_index(int nobs,double q)
{
  // Given nobs in order, find index 0:(nobs-1) with quantile q
  int qi;
  qi=(int) ceil(((double) nobs)*q);
  
  
  return(qi);
}

double get_summary_basic(int i,int j,int m,double q)
{
  double r;
  int qi;
  
  switch(daop.type){

  case 1:
    r=daop.hs[j].s[m][i];
    break;

  case 2:
    if(daop.hs[j].xindx==(-1))
      {
	// Split on number of obs
	r=daop.hs[j].s[m][i];
      }else{
      // split on quantile 
      if(daop.hs[j].nobs[m][i]==0){
	r=0;
      }else{
	// find quantile (note: s[] is sorted in descending order (smallest last))
	qi=quantile_index(daop.hs[j].nobs[m][i],q)-1;
	r=daop.hs[j].s[m][daop.hs[j].start_vec[m][i]+daop.hs[j].nobs[m][i]-1-qi];

      }
    }
    break;

  case 3:   // split on order-statistic  
        if(daop.hs[j].xindx==(-1))
      {
	// Split on number of obs
	r=daop.hs[j].s[m][i];
      }else{
      // split on quantile 
      if(daop.hs[j].nobs[m][i]==0){
	r=0;
      }else{
	// extract order statistic, if does not exist, return a small number
	qi=(int) q;
	if(qi>daop.hs[j].nobs[m][i])
	  {
	    r=-BIG_NUMBER;
	  }else{
	  qi=qi-1;
	  r=daop.hs[j].s[m][daop.hs[j].start_vec[m][i]+daop.hs[j].nobs[m][i]-1-qi];
	}

      }
    }

  }
  return(r);
}

double get_summary_basic2(int i,int j,int m,int m2,double q)
{
  // Modification for windowed summaries 
  
  double r;
  int qi,num_obs;
  int type_mod,k,nobs,counter;

   
  type_mod=daop.type;

  if(type_mod==5&&m2<0)
    type_mod=2;

  if(type_mod==6&&m2<0)
    type_mod=3;

  if(type_mod==4&&m2<0)
    type_mod=1;
  
  
 
  switch(type_mod){
  
  case 1:
    r=daop.hs[j].s[m][i];
    break;

  case 2:
    if(daop.hs[j].xindx==(-1))
      {
	// Split on number of obs
	r=daop.hs[j].s[m][i];
      }else{
      // split on quantile 
      if(daop.hs[j].nobs[m][i]==0){
	r=0;
      }else{
	// find quantile (note: s[] is sorted in descending order (smallest last))
	qi=quantile_index(daop.hs[j].nobs[m][i],q)-1;
	r=daop.hs[j].s[m][daop.hs[j].start_vec[m][i]+daop.hs[j].nobs[m][i]-1-qi];

      }
    }
    break;

  case 3:   // split on order-statistic  
        if(daop.hs[j].xindx==(-1))
      {
	// Split on number of obs
	r=daop.hs[j].s[m][i];
      }else{
      // split on quantile 
      if(daop.hs[j].nobs[m][i]==0){
	r=0;
      }else{
	// extract order statistic, if does not exist, return a small number
	qi=(int) q;
	if(qi>daop.hs[j].nobs[m][i])
	  {
	    r=-BIG_NUMBER;
	  }else{
	  qi=qi-1;
	  r=daop.hs[j].s[m][daop.hs[j].start_vec[m][i]+daop.hs[j].nobs[m][i]-1-qi];
	}

      }
	}
	break;

  case 4: // Window summary ...
    // get number of observations
    if(m2<0)
      {
	// ... ignore m2.
	r=daop.hs[j].s[m][i];
      }else{
      num_obs=daop.hs[daop.nvar_history].s[m][i]-daop.hs[daop.nvar_history].s[m2][i];
      if(j==daop.nvar_history)
	{
	  r=num_obs;
	}else{
      if(num_obs>.001) // number of obs between m and m2
	{
	  r=(daop.hs[j].s[m][i]*daop.hs[daop.nvar_history].s[m][i]-daop.hs[j].s[m2][i]*daop.hs[daop.nvar_history].s[m2][i])/num_obs;
	}else{
	r=0;
      }
      }
    }
       
    break;


  case 5: // Window quantile summary
    if(daop.hs[j].xindx==(-1))
      {
      
	// Split on number of obs
	//	r=(double) (daop.hs[j].nobs[m][i]-daop.hs[j].nobs[m2][i]);
		r=(double) (daop.hs[j].s[m][i]-daop.hs[j].s[m2][i]);

      }else{
      // split on quantile
      nobs=(daop.hs[j].nobs[m][i]-daop.hs[j].nobs[m2][i]);
 
      if(nobs<=0){
	r=0;
      }else{
	// find quantile (note: s[] is sorted in descending order (smallest last))
	qi=quantile_index(nobs,q)-1;
	// r=daop.hs[j].s[m][daop.hs[j].start_vec[m][i]+daop.hs[j].nobs[m][i]-1-qi];
	counter=0;
	for(k=0;k<daop.hs[j].nobs[m][i];k++)
	  {
	    if(daop.hs[j].deltav[m][daop.hs[j].start_vec[m][i]+daop.hs[j].nobs[m][i]-1-k]>daop.delta[m2])
	      counter++;

		if(counter==(qi+1))
		{
		  r=daop.hs[j].s[m][daop.hs[j].start_vec[m][i]+daop.hs[j].nobs[m][i]-1-k];
		  break;
		}
	  }
	
	if(counter==0)
	  r=0;
	
	// 
	
      }
    }
    break;


  case 6: // Window order stat 
    if(daop.hs[j].xindx==(-1))
      {
      
	// Split on number of obs
	//	r=(double) (daop.hs[j].nobs[m][i]-daop.hs[j].nobs[m2][i]);
		r=(double) (daop.hs[j].s[m][i]-daop.hs[j].s[m2][i]);

      }else{
      // split on orderstat
      nobs=(daop.hs[j].nobs[m][i]-daop.hs[j].nobs[m2][i]);

	qi=(int) q;
	if(qi>nobs)
	  {
	    r=-BIG_NUMBER;
	  }else{
	

	counter=0;
	for(k=0;k<daop.hs[j].nobs[m][i];k++)
	  {
	    if(daop.hs[j].deltav[m][daop.hs[j].start_vec[m][i]+daop.hs[j].nobs[m][i]-1-k]>daop.delta[m2])
	      counter++;

		if(counter==(qi+1))
		{
		  r=daop.hs[j].s[m][daop.hs[j].start_vec[m][i]+daop.hs[j].nobs[m][i]-1-k];
		  break;
		}
	  }
	if(counter==0)
	  r=0;
	  
	}
      
 
    }
    break; 


    
  }
  return(r);
}


int sample_m2(int m)
{
  int m2;
  double r,md;
  //GetRNGstate();  
  if(m==0)
    {
      m2=-99;
    }else{
    md=(double) m;
    r=unif_rand();
    m2=(int) floor(md*r);
    if(m2==0)
      m2=-99;
  }
  //PutRNGstate();
  return(m2);
}


int qtry_get(int j,int m,int m2)
{
  // Get number of quantiles/order-statistics for splitting
  // .. also, sample these if possible number exceeds daop.qtry
  
  
  int qtry_use;
  

  	  if(j<daop.nvar_history){
	    if(daop.type==2||(daop.type==5&&m2<0))
	    qtry_use=sample_quantiles2(daop.quantile_sample,daop.hs[j].max_nobs[m]);
	  // should sample quantiles/order statistics here ... (number of possible will depend on delta/delta0)

	  //Rprintf("sampling order-stat maxobs=%d \n",daop.hs[j].max_nobs[m]);
	    if(daop.type==3||(daop.type==6&&m2<0))
	    qtry_use=sample_orderstat(daop.order_sample,daop.hs[j].max_nobs[m]);

	   if(daop.type==4||daop.type==1)
	     qtry_use=1;


	   if(daop.type==5&&m2>(-1)) // windowed quantile 
	     qtry_use=sample_quantiles_w(daop.quantile_sample,m,m2); 


	   if(daop.type==6&&(m2>(-1))) // windowed quantile 
	     qtry_use=sample_orderstat_w(daop.order_sample,m,m2); 

	  // should sample quantiles/order statistics here ... (number of possible will depend on delta/delta0)
	  }else{
	    qtry_use=1;
	  }

	  


	  return(qtry_use);

}


void split_history_basic(int nodeid,struct split_result *best_split)
{
  // split node based on historical summaries


 
  
  struct split_result sr;
  double *x,*y;
  int *xi;
  double q;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j,l,m,u;
  // int predictors[1000];
  int qtry_use,mtry_use;
  int m2;
  int *predictors;

  predictors=(int*) malloc(sizeof(int)*daop.n*2);

    
  x=daop.daux5;
  y=daop.daux4;
  xi=daop.iaux5;

  init=0;
  best_split->success=0; // until proven otherwise....

  //if(daop.mtry<daop.nvar_history){ //
  // sample predictors 
  counter=0;
  for(k=0;k<=daop.nvar_history;k++)
    {
 
	  predictors[counter]=k;
	  counter++;
	
    }

  permute(predictors,(counter));
  //}

  mtry_use=daop.nvar_history+1;
  if(daop.mtry<mtry_use)
    mtry_use=daop.mtry;
  
  // Sample delta 
  sample_delta_index(daop.delta_sample);


    
  if(daop.type==1||daop.type==4)
    daop.qtry=1;
  
  
  //for(j=0;j<=daop.nvar_history;j++)
  for(u=0;u<mtry_use;u++)
    {
      j=predictors[u];
      
       qtry_use=daop.qtry;
      m2=-99;
      for(k=0;k<daop.dtry;k++)
	{

	  m=daop.delta_sample[k];

	  // ... sample number 0:(m-1)
	 
	  if(daop.type>3)
	     m2=sample_m2(m);


	  qtry_use=qtry_get(j,m,m2);


	  for(l=0;l<qtry_use;l++)
	    {

	      if(j<daop.nvar_history){
	      if(daop.type==2||daop.type==5)
		q=daop.quantile_sample[l];


	    
	      if(daop.type==1||daop.type==4)
		q=0;
	      
	      if(daop.type==3||daop.type==6)
		q=(double) daop.order_sample[l];
	      }
	  // get data
	  HASH_FIND_INT(nodes,&nodeid,s);

	  counter=0;
	  if(s==NULL){}else{ 

	    HASH_ITER(hh,s->rows,r,rtemp)
	      {
		i=r->row_number;
		if(daop.train[i]==1){
		  //y[counter]=daop.x[daop.yindx][i];
		  y[counter]=daop.response[i];
		  x[counter]=get_summary_basic2(i,j,m,m2,q);
		  counter++;
		}
	      }
	  }
	  
	  if(counter>0)
	    {
	      if(daop.ncat<0){
		basic_split(x,y,counter,xi,daop.min_nodesize,&sr);
	      }else{
		basic_split_gini(x,y,counter,xi,daop.min_nodesize,&sr);
	      }
	      
	      if(sr.success==1)
		{
		  if(init==0||((best_split->error)>sr.error))
		    {
		      best_split[0]=sr;
		      best_split->varindx=daop.hs[j].xindx;  // 
		      best_split->delta=daop.delta[m];//daop.hs[j].delta[m];
		      best_split->delta0=daop.delta0[m];//daop.hs[j].delta0[m];
		      best_split->scolumn=(double) m;
		      best_split->scolumn2=(double) m2;
		      best_split->hsindx=(double) j;
		      if((daop.type==2||daop.type==3||daop.type==5||daop.type==6)&&(daop.hs[j].xindx!=-1))
			best_split->tau=q;

		      if(daop.type>3&&m2>(-1))
			best_split->delta0=daop.delta[m2]+.0001;//daop.hs[j].delta[m2]+.0001; // SOME SMALL NUMBER ... beware 

		      
		      init=1;
		    }
		  
		}
	      

	    }
	    }
	}
    }
  free(predictors);
}





void print_concurrent()
{
  int k;
  for(k=0;k<daop.nvar_concurrent;k++)
    Rprintf("k=%d: %d ",k,daop.splitvar_concurrent[k]);
  Rprintf("\n");

}




// -------------------------------------------------------------------
// ---------------- BASIC  -------------------------------------------
// -------------------------------------------------------------------




void update_nodevector_basic(int nodeid,int nodeL,int nodeR,struct split_result *sr)
{
  //void tsummary(double *nrec,double *nrec_condition,double cut,double delta,int vindx)
  struct data_options *d;
  int i;
  double ff;
  struct node *s;
  struct rnumber *r,*rtemp;
  double *x;
  double cut;
  int hsplit,hsindx,scolumn,scolumn2;
  double q;
  
    d=&daop;

    hsplit=0;
  if(sr->delta<0) // c-split
    {
      x=d->x[sr->varindx];
    }else{ // h-split
    hsplit=1;
    hsindx=sr->hsindx;
    scolumn=sr->scolumn;
    q=sr->tau;
    scolumn2=sr->scolumn2; 
  }

  cut=sr->point;
  
  HASH_FIND_INT(nodes,&nodeid,s);
  if(s==NULL){
    Rprintf(" Updating node list, node %d does not exist ...\n",nodeid);
  }else{

    // for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {
      i=r->row_number;
      //ff=d->x[vindx][i];
      if(hsplit==0){
      ff=x[i];
      }else{
	ff=get_summary_basic2(i,hsindx,scolumn,scolumn2,q);
	
	}
	  if(ff<cut)
	    {
	      // left node 
	      add_row(nodeL,i);
	      //d->node[i]=nodeL;
	    }else{
	    // right node 
	    add_row(nodeR,i);
	    //d->node[i]=nodeR;
	  }
	
    }
  

    delete_node(nodeid);
  }
}






void update_treematrix_basic(struct split_result *sr,int nodeid,int node_row,int node_counter)
{


  int m,ni,k;
  double a;
  double avec[MAXCAT];
  struct data_options *d;
  d=&daop;

  
  // update split info for nodeid (node just split)
  d->tree_matrix[4][node_row]=((double) sr->varindx); // split variable
  d->tree_matrix[5][node_row]=sr->tau; // fraction where split is performed 
  d->tree_matrix[6][node_row]=sr->point; //  variable cut
  d->tree_matrix[7][node_row]=sr->delta; // delta
  d->tree_matrix[8][node_row]=(d->row_counter+1); // row number corresponding left node 
  d->tree_matrix[9][node_row]=(d->row_counter+2); // row number corresponding right node
  d->tree_matrix[10][node_row]=sr->error_reduction;
  d->tree_matrix[11][node_row]=sr->nL+sr->nR;
  d->tree_matrix[12][node_row]=sr->delta0;
  d->tree_matrix[13][node_row]=sr->hsindx;
  d->tree_matrix[14][node_row]=sr->scolumn;
  d->tree_matrix[15][node_row]=sr->scolumn2;

  


						     
  // left node prediction
  ni=node_counter-1;
  a=node_prediction_uthash(ni);
  a=d->lambda*a;
  //Rprintf(" prediction left node=%lf \n",a);
  d->row_counter++;

  // add to tree matrix 
  d->tree_matrix[0][d->row_counter]=((double) d->tree_counter); // tree indicator
  d->tree_matrix[1][d->row_counter]=((double) nodeid); // parent node id
  d->tree_matrix[2][d->row_counter]=((double) (node_counter-1)); // node id
  d->tree_matrix[3][d->row_counter]=a; // prediction
  for(m=4;m<=15;m++)
    d->tree_matrix[m][d->row_counter]=-99; // split info
  d->tree_matrix[11][d->row_counter]=sr->nL; 

  // add class probabilities 
 if(daop.ncat>0)
   {
     node_prediction_gini(ni,avec);
     for(k=0;k<daop.ncat;k++)
       d->tree_matrix[16+k][d->row_counter]=avec[k];

   }
  

  // right node prediction
  a=node_prediction_uthash((node_counter));
  a=d->lambda*a;
  d->row_counter++;
  // add to tree matrix 
  d->tree_matrix[0][d->row_counter]=((double) d->tree_counter); // tree indicator
  d->tree_matrix[1][d->row_counter]=((double) nodeid); // parent node id
  d->tree_matrix[2][d->row_counter]=((double) (node_counter)); // node id
  d->tree_matrix[3][d->row_counter]=a; // prediction
  for(m=4;m<=15;m++)
    d->tree_matrix[m][d->row_counter]=-99; // split info 
  d->tree_matrix[11][d->row_counter]=sr->nR; // prediction


  // add class probabilities 
 if(daop.ncat>0)
   {
     
     node_prediction_gini(node_counter,avec);
     for(k=0;k<daop.ncat;k++)
        d->tree_matrix[16+k][d->row_counter]=avec[k];

   }

}






void tree(int *ptr_nsplit)
{

  int i,m,k,j;
  int nsplit,node_counter;
  int selector;
  int best_var,best_var_2;
  int success,success_2,nodeid;
  int *nodes_to_split;
  int *node_rownumber;
  int nodes_not_split;
  double a;
  struct data_options *d;
  int ni;
  struct split_result best_split;
  struct split_result best_split_time;
  struct split_result sr;
  double ssq_reduction,ssq_reduction_time;
  
  d=&daop;
  
  nsplit=ptr_nsplit[0];
  nodes_to_split=(int*) malloc(sizeof(int)*10*nsplit);
  node_rownumber=(int*) malloc(sizeof(int)*10*nsplit);
  
  d->row_counter++; // increment  
  d->tree_start[d->tree_counter]=d->row_counter;

  // initialize root node 
  j=0;
  for(i=0;i<d->n;i++)
    add_row(j,i);

  // Set up root-node info 
  nodes_to_split[0]=0;
  node_rownumber[0]=d->row_counter;
  // mean of response
  ni=0;
  a=node_prediction_uthash(ni); 
  a=d->lambda*a;       // shrink mean by shrinkage factor lambda.

  // tree matrix entry for root node
  d->tree_matrix[0][d->row_counter]=d->tree_counter; // tree indicator
  d->tree_matrix[1][d->row_counter]=-99; // parent node id
  d->tree_matrix[2][d->row_counter]=0; // node id
  d->tree_matrix[3][d->row_counter]=a; // prediction
  d->tree_matrix[4][d->row_counter]=-99; // .. this is the terminal node for now

 
 
  
  nodes_not_split=1;
  node_counter=0;
  nodeid=0;
  for(k=0;k<nsplit;k++)
    {

      if(k>node_counter)
	break;
      
      nodeid=nodes_to_split[k];

 
      // split node
      split_standard(nodeid,&best_split);
      success=best_split.success;
      if(success==1){
	error_reduction_split_result(&best_split);
      }

 
      if(daop.nvar_history>0){
      split_history_basic(nodeid,&best_split_time);

      if(best_split_time.success)
	{
	  error_reduction_split_result(&best_split_time);
	  if((success==1&&(best_split_time.error_reduction>best_split.error_reduction))||(success==0)){
	    best_split=best_split_time;
	    success=1;
	  }
	}
      }
      
      if(success==1)
	{
	  for(i=1;i<=2;i++){
	    node_counter++;
	    nodes_to_split[node_counter]=node_counter;
	    node_rownumber[node_counter]=(d->row_counter+i);
	  }


	  update_nodevector_basic(nodeid,(node_counter-1),node_counter,&best_split);
	  update_treematrix_basic(&best_split,nodeid,node_rownumber[nodeid],node_counter); 


	}


 
    }

  d->tree_end[d->tree_counter]=d->row_counter;
  free(nodes_to_split); 
  free(node_rownumber); 

}







void tree_gini(int *ptr_nsplit)
{

  int i,m,k,j;
  int nsplit,node_counter;
  int selector;
  int best_var,best_var_2;
  int success,success_2,nodeid;
  int *nodes_to_split;
  int *node_rownumber;
  int nodes_not_split;
  double a;
  struct data_options *d;
  int ni;
  struct split_result best_split;
  struct split_result best_split_time;
  struct split_result sr;
  double ssq_reduction,ssq_reduction_time;
  
  d=&daop;
  
  nsplit=ptr_nsplit[0];
  nodes_to_split=(int*) malloc(sizeof(int)*10*nsplit);
  node_rownumber=(int*) malloc(sizeof(int)*10*nsplit);
  
  d->row_counter++; // increment  
  d->tree_start[d->tree_counter]=d->row_counter;

  // initialize root node 
  j=0;
  for(i=0;i<d->n;i++)
    add_row(j,i);

  // Set up root-node info 
  nodes_to_split[0]=0;
  node_rownumber[0]=d->row_counter;
  // mean of response
  ni=0;
  a=node_prediction_uthash(ni); 
  a=d->lambda*a;       // shrink mean by shrinkage factor lambda.

  // tree matrix entry for root node
  d->tree_matrix[0][d->row_counter]=d->tree_counter; // tree indicator
  d->tree_matrix[1][d->row_counter]=-99; // parent node id
  d->tree_matrix[2][d->row_counter]=0; // node id
  d->tree_matrix[3][d->row_counter]=a; // prediction
  d->tree_matrix[4][d->row_counter]=-99; // .. this is the terminal node for now

 
 
  
  nodes_not_split=1;
  node_counter=0;
  nodeid=0;
  for(k=0;k<nsplit;k++)
    {
     
      if(k>node_counter)
	break;
      
      nodeid=nodes_to_split[k];

  
      // split node
      split_standard_gini(nodeid,&best_split);
     success=best_split.success;
      if(success==1){
	error_reduction_split_gini(&best_split);
      }


      if(daop.nvar_history>0){
      split_history_basic(nodeid,&best_split_time);
 
      if(best_split_time.success)
	{
	  error_reduction_split_gini(&best_split_time);
	  if((success==1&&(best_split_time.error_reduction>best_split.error_reduction))||(success==0)){
	    best_split=best_split_time;
	    success=1;
	  }
	}
      }
      
      if(success==1)
	{
	  for(i=1;i<=2;i++){
	    node_counter++;
	    nodes_to_split[node_counter]=node_counter;
	    node_rownumber[node_counter]=(d->row_counter+i);
	  }


	  update_nodevector_basic(nodeid,(node_counter-1),node_counter,&best_split);
	  update_treematrix_basic(&best_split,nodeid,node_rownumber[nodeid],node_counter); 


	}


 
    }
 
  d->tree_end[d->tree_counter]=d->row_counter;
  free(nodes_to_split);
  free(node_rownumber);

}






void read_basic(double *x,int *id,double *time,int *n,int *p,int *yindx,int *nsplit,int *nboost,int *mtry,int *nodesize,
		 int *concurrent,int *nconcurrent,  
		int *historic,int *nhistoric)
{
  int i,j,ncol,k,counter;
  struct data_options *d;
  double *v;
  int max_history_length;
  int nsampw;
  d=&daop;

  daop.s_allocated=0;
  daop.method=1;
  
  daop.ncc=0;

  daop.ncat=-10; // not classification, until set. 
 
  daop.yindx=yindx[0];
  daop.nsplit=n[0];
  
  max_history_length=1;
  // -------------------------------------------------------
  ncol=15+MAXCAT;  // incremented for window-summaries (prior was 14)


 
  
  

  d->mtry=mtry[0];
  d->lambda=1; // shrinkage parameter 
  //  d->nsplit=nsplit[0];
  d->nboost=nboost[0];
  // tree matrix
 

  
  d->row_counter=-1; // this is incremented as rows are added
  d->tree_counter=0;

  d->tree_start=(int*) malloc(sizeof(int)*d->nboost); //fr
  d->tree_end=(int*) malloc(sizeof(int)*d->nboost); //fr

      d->subsample=1;
      d->sample_fraction=.5;

 

  d->min_nodesize=nodesize[0];
  
  d->p=p[0];
  d->n=n[0];

  // Should not reallocate, just copy addresses.. 
  d->x=(double**) malloc(sizeof(double*)*d->p);  //fr
  for(i=0;i<d->p;i++)
      d->x[i]=&(x[(d->n)*i]);



  d->nvar_history=nhistoric[0];
  d->splitvar_history=historic;
  d->nvar_concurrent=nconcurrent[0];
  d->splitvar_concurrent=concurrent;

  if(d->nvar_history)
    d->hs=(struct histvar*) malloc(sizeof(struct histvar)*(d->nvar_history+1));
 


  
  // response 
  d->response=(double*) malloc(sizeof(double)*d->n);
  for(i=0;i<d->n;i++)
    d->response[i]=d->x[yindx[0]][i];  //fr


  d->id=(int*) malloc(sizeof(int)*d->n);
  for(i=0;i<d->n;i++)
    d->id[i]=id[i];

  d->time=time;

  
  d->train=(int*) malloc(sizeof(int)*d->n);//fr
     for(j=0;j<d->n;j++)
    d->train[j]=1;

   
  
  d->predictions=(double*) malloc(sizeof(double)*d->n);//fr
  for(j=0;j<d->n;j++)
    d->predictions[j]=0;


// oob counter vector 
    d->oob=(int*) malloc(sizeof(int)*d->n);//fr
  for(j=0;j<d->n;j++)
    d->oob[j]=0;
  
  // node vector
  d->node=(int*) malloc(sizeof(int)*d->n);//fr
  for(i=0;i<d->n;i++)
    d->node[i]=0;

  
  // scratch space  fr-all
  d->daux1=(double*) malloc(sizeof(double)*d->n*max_history_length);
  d->daux2=(double*) malloc(sizeof(double)*d->n*max_history_length);
  d->daux3=(double*) malloc(sizeof(double)*d->n*max_history_length);
  d->daux4=(double*) malloc(sizeof(double)*d->n*max_history_length);
  d->daux5=(double*) malloc(sizeof(double)*d->n*max_history_length);
  d->daux6=(double*) malloc(sizeof(double)*d->n*max_history_length);

  d->iaux1=(int*) malloc(sizeof(int)*d->n*max_history_length);
  d->iaux2=(int*) malloc(sizeof(int)*d->n*max_history_length);
  d->iaux3=(int*) malloc(sizeof(int)*d->n*max_history_length);
  d->iaux4=(int*) malloc(sizeof(int)*d->n*max_history_length);
  d->iaux5=(int*) malloc(sizeof(int)*d->n*max_history_length);
  d->iaux6=(int*) malloc(sizeof(int)*d->n*max_history_length);
  d->iaux7=(int*) malloc(sizeof(int)*d->n*max_history_length);


  d->LARGE_NUMBER=999999;

 
}


 

void update_residual_basic(int tree_indx)
{
  int k,i,start,go,ni,ti,end;
  int ntnodes;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  int *tnodes;
  double *node_predictions;

  tnodes=(int*) malloc(sizeof(int)*daop.n*2);
  node_predictions=(double*) malloc(sizeof(double)*daop.n*2);

  
  ntnodes=0;
  start=daop.tree_start[tree_indx];
  end=daop.tree_end[tree_indx];
  
  for(k=start;k<=end;k++)
    {

      ti=(int) daop.tree_matrix[0][k];
      if(ti==tree_indx)
	{
	  ni=((int) daop.tree_matrix[2][k]);
	  //Rprintf(" node id =%d k=%d \n",ni,k);
	  node_predictions[ni]=daop.tree_matrix[3][k];
	  if(daop.tree_matrix[4][k]<(-1))
	    {
	      tnodes[ntnodes]=ni;
	      ntnodes++;
	    }
	}else{
	go=0;
	break;
      }

    }

  for(k=0;k<ntnodes;k++)
    {

      ni=tnodes[k];
      HASH_FIND_INT(nodes,&ni,s);
      if(s==NULL){}else{

	HASH_ITER(hh,s->rows,r,rtemp){
      
	  i=r->row_number;
      if(daop.rf!=1)
	{
	  daop.response[i]=(daop.response[i]-node_predictions[ni]);
	  daop.predictions[i]+=(node_predictions[ni]);
	}else{

	if(daop.train[i]==0) // out-of-bag
	  {
	    daop.predictions[i]+=(node_predictions[ni]);
	    daop.oob[i]++;
	  }
      }

      }
      }
    }


  free(tnodes);
  free(node_predictions);

}





void update_residual_gini(int tree_indx,int *oob)
{
  int k,i,start,go,ni,ti,end,m;
   int ntnodes;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  int NCOL_TREE;
  double *aa;
  int *tnodes;
  double **node_predictions;
  
  NCOL_TREE=16;

  node_predictions=(double**) malloc(sizeof(double*)*daop.n);
  aa=(double*) malloc(sizeof(double)*daop.n*daop.ncat);
  for(k=0;k<daop.n;k++)
    node_predictions[k]=&(aa[k*daop.ncat]);
  tnodes=(int*) malloc(sizeof(int)*daop.n);

  
  ntnodes=0;
  start=daop.tree_start[tree_indx];
  end=daop.tree_end[tree_indx];
  
  for(k=start;k<=end;k++)
    {

      ti=(int) daop.tree_matrix[0][k];
      if(ti==tree_indx)
	{
	  ni=((int) daop.tree_matrix[2][k]);
	  for(m=0;m<daop.ncat;m++)
	    node_predictions[ni][m]=daop.tree_matrix[NCOL_TREE+m][k];
	  
	  if(daop.tree_matrix[4][k]<(-1))
	    {
	      tnodes[ntnodes]=ni;
	      ntnodes++;
	    }
	}else{
	go=0;
	break;
      }

    }

  for(k=0;k<ntnodes;k++)
    {

      ni=tnodes[k];
      HASH_FIND_INT(nodes,&ni,s);
      if(s==NULL){}else{

	HASH_ITER(hh,s->rows,r,rtemp){
      
	  i=r->row_number;

	  oob[i]=0;
	if(daop.train[i]==0) // out-of-bag
	  {
	    for(m=0;m<daop.ncat;m++)
	      daop.predictions_gini[i][m]+=(node_predictions[ni][m]);
	    daop.oob[i]++;

	    oob[i]=1;
	  }
     

      }
      }
    }



  free(node_predictions[0]);
  free(node_predictions);
  free(tnodes);

}




void update_residual_basic2(int tree_indx,int *oob,double *pred)
{
  // This function extracts oob vector and tree predictions for
  // this tree, useful in parallel runs.
  
  int k,i,start,go,ni,ti,end;
  // double node_predictions[10000]; // .. beware.. 
  //int tnodes[10000];// terminal nodes
  double pr;
  int ntnodes;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  double *node_predictions;
  int *tnodes;

  tnodes=(int*) malloc(sizeof(int)*2*daop.n);
  node_predictions=(double*) malloc(sizeof(double)*2*daop.n);
  

  ntnodes=0;
  start=daop.tree_start[tree_indx];
  end=daop.tree_end[tree_indx];
  
  for(k=start;k<=end;k++)
    {

      ti=(int) daop.tree_matrix[0][k];
      if(ti==tree_indx)
	{
	  ni=((int) daop.tree_matrix[2][k]);
	  node_predictions[ni]=daop.tree_matrix[3][k];
	  if(daop.tree_matrix[4][k]<(-1))
	    {
	      tnodes[ntnodes]=ni;
	      ntnodes++;
	    }
	}else{
	go=0;
	break;
      }

    }

  for(k=0;k<ntnodes;k++)
    {

      ni=tnodes[k];
      HASH_FIND_INT(nodes,&ni,s);
      if(s==NULL){}else{

	HASH_ITER(hh,s->rows,r,rtemp){
      
	  i=r->row_number;
      if(daop.rf!=1)
	{
	  if(daop.family==1){
	    daop.response[i]=(daop.response[i]-node_predictions[ni]);
	    daop.predictions[i]+=(node_predictions[ni]);
	  }
	  if(daop.family==2){
	    // daop.response[i]=(daop.response[i]-node_predictions[ni]);
	    daop.predictions[i]+=(node_predictions[ni]);
	    pr=1/(1+exp(-daop.predictions[i]));
	    daop.response[i]=(daop.x[daop.yindx][i]-pr);


	  }
	  
	}else{

	pred[i]=node_predictions[ni];
	oob[i]=0;
	if(daop.train[i]==0) // out-of-bag
	  {
	    daop.predictions[i]+=(node_predictions[ni]);
	    daop.oob[i]++;
	    oob[i]=1;
	  }
      }

      }
      }
    }


  free(tnodes);
  free(node_predictions);
  

}








void free_basic()
{
  int i,j,k;
  struct data_options *d;

  d=&daop;

   delete_oob();

 free(d->predictions);
 
 free(d->train); 
   free(d->tree_start);
  free(d->tree_end);

  free(d->x);


  free(d->response);


  free(d->id);
  free(d->oob);

  free(d->node);
    // scratch space
  free(d->daux1); 
  free(d->daux2); 
  free(d->daux3); 
  free(d->daux4); 
  free(d->daux5); 
  free(d->daux6);

 
  free(d->iaux1);
  free(d->iaux2);
  free(d->iaux3);
  free(d->iaux4);
  free(d->iaux5);
  free(d->iaux6);
  free(d->iaux7);


  if(daop.nvar_history>0)
    {
      free(daop.delta_sample);

      if(daop.s_allocated==1)
	{
	  for(k=0;k<=daop.nvar_history;k++)
	    {
	      for(j=0;j<daop.ndelta;j++){
		free(daop.hs[k].s[j]);
		//
		if((daop.type==2||daop.type==3||daop.type==5||daop.type==6)&&(k<daop.nvar_history))
		  {
		    free(daop.hs[k].nobs[j]);
		    free(daop.hs[k].start_vec[j]);
		    free(daop.hs[k].deltav[j]);
		    
		  }
	      }

	      free(daop.hs[k].s);
	      //
	      if((daop.type==2||daop.type==3||daop.type==5||daop.type==6)&&(k<daop.nvar_history)){
		free(daop.hs[k].start_vec);
		free(daop.hs[k].nobs);
		free(daop.hs[k].max_nobs);
		free(daop.hs[k].deltav);
		
	      }
	    }

	}
      if(daop.type==3||daop.type==6)
	free(daop.order_sample);

      if(daop.type==2||daop.type==5)
	free(daop.quantile_sample);
      
      free(daop.hs);
    }
  
  free(daop.tree_matrix);
}

void mtree(int *ptr_nboost,int *nsplit,double *te,double *predictions,double *trees)
{
  int b,i,j;
  int B;
  int NCOL_TREE;

 
  NCOL_TREE=16;

  if(daop.ncat>0){ // 
    NCOL_TREE=NCOL_TREE+daop.ncat;

    daop.predictions_gini=(double**) malloc(sizeof(double*)*daop.n);
    for(j=0;j<daop.n;j++)
      daop.predictions_gini[j]=&(predictions[daop.ncat*j]);
    
  }
    
  GetRNGstate();
  B=ptr_nboost[0];

  // DEBUG

  daop.method=1;
  //daop.rf=1;
  daop.subsample=1;
  daop.row_counter=-1;
  daop.nsplit=nsplit[0];
  // .. only active for logistic regression
  //  initialize_prediction_response();



  // ... allocate tree matrix
  daop.tree_matrix=(double**) malloc(sizeof(double*)*NCOL_TREE);
  for(i=0;i<NCOL_TREE;i++)
    daop.tree_matrix[i]=&(trees[B*(daop.nsplit*2+1)*i]);

  if(daop.rf==1){
    delete_oob();   
    set_id_rows(daop.id);
  }

 
  // Ensemble loop   
  for(b=0;b<B;b++)
    {

      
      if(daop.subsample==1&&daop.rf==1)
	subsample_id(b);
       

      if(daop.ncat<0){
	tree(&(daop.nsplit));
      }else{
	tree_gini(&(daop.nsplit));
      }
 
      update_terminal_nodes_wrapper(daop.tree_counter);
      
      if(daop.ncat<0){
	update_residual_basic2(daop.tree_counter,
			       &(daop.oob_tree[daop.n*b]),&(daop.pred_tree[daop.n*b]));
      }else{
	update_residual_gini(daop.tree_counter,&(daop.oob_tree[daop.n*b]));
      }

      // update_residual_wrapper(daop.tree_counter);
      daop.tree_counter++;

      if(daop.rf==1){
	if(daop.ncat<0){
	  te[b]=test_error();
	}else{
	  te[b]=test_error_gini();
	  Rprintf(" error=%lf ",te[b]);
	}
      }

      if(daop.rf!=1)
	te[b]=train_error(); //test_error_wrapper();
	
      
      delete_all_nodes();      

    }

  
   
 
  //free(daop.tree_matrix);
  
  //Rprintf("Total loop count=%d \n",daop.ncc);
  if(daop.ncat<0){
  for(i=0;i<daop.n;i++)
   predictions[i]=daop.predictions[i];
  }

  if(daop.ncat)
    free(daop.predictions_gini);
  
   PutRNGstate();
}









void ctree(double *x,int *id,double *time,int *n,int *p,int *yindx,int *nsplit,int *nboost,int *mtry,int *nodesize,
	   int *concurrent,int *nconcurrent,  
	   int *historic,int *nhistoric,
	   double *delta,double *delta0,int *ndelta,int *type,double *quantiles,int *nquantiles,int *dtry,int *qtry,int *rf,double *lambda,int *family,int *ncat,double *sample_fraction,
	   double *te,double *predictions,double *trees,int *oob,double *tpred,int *nrows
	   )
{


  daop.qtry=qtry[0];
  if(type[0]==2||type[0]==5)
    daop.quantile_sample=(double*) malloc(sizeof(double)*nquantiles[0]); // changed from nquantiles[0]

  daop.family=family[0];
  daop.oob_tree=oob;
  daop.pred_tree=tpred;
  
  read_basic(x,id,time,n,p,yindx,nsplit,nboost,mtry,nodesize,concurrent,nconcurrent,historic,nhistoric);
  create_history(delta,delta0,ndelta,type,quantiles,nquantiles,dtry);
  daop.ncat=ncat[0]; 

  daop.sample_fraction=sample_fraction[0];
  // Set rf or not
  daop.rf=rf[0];
  daop.lambda=lambda[0];
     
  mtree(nboost,nsplit,te,predictions,trees);
  nrows[0]=daop.row_counter+1;
   
  free_basic();
 

}


void ctree_predict(double *x,int *id,double *time,int *n,int *p,int *yindx,int *nboost,
	   int *concurrent,int *nconcurrent,  
	   int *historic,int *nhistoric,
		   double *delta,double *delta0,int *ndelta,int *type,int *ncat,double *quantiles,int *nquantiles,
		   double *trees,int *nrows,int *all_trees,double *predictions
	   )
{
  int NCOL_TREE;
  int mtry,nsplit,nodesize,dtry;
  int k,curr_tree;
  int num_trees,counter,offset;

  
  NCOL_TREE=16;
  mtry=1;
  nsplit=1;
  nodesize=1;
  dtry=1;


  daop.qtry=1;
  if(type[0]==2||type[0]==5)
    daop.quantile_sample=(double*) malloc(sizeof(double)*nquantiles[0]);

 
 
    read_basic(x,id,time,n,p,yindx,&nsplit,nboost,&mtry,&nodesize,concurrent,nconcurrent,historic,nhistoric);

  free(daop.tree_start);

   create_history(delta,delta0,ndelta,type,quantiles,nquantiles,&dtry);

 
    daop.ncat=ncat[0];
    if(daop.ncat>0)
      NCOL_TREE=NCOL_TREE+daop.ncat;
 
    // .. trees
    daop.tree_matrix=(double**) malloc(sizeof(double*)*NCOL_TREE);
    for(k=0;k<NCOL_TREE;k++)
      daop.tree_matrix[k]=&(trees[k*nrows[0]]);

 

    
  
      // prev   if(0){   
     curr_tree=0;
     //    daop.tree_start[0]=curr_tree;
    // count number of trees
    num_trees=1;
      for(k=0;k<nrows[0];k++)
      {
	if(daop.tree_matrix[0][k]!=curr_tree)
	  {
	    curr_tree=daop.tree_matrix[0][k];
	    num_trees++;
	  }
      }
   
     daop.nboost=num_trees;
    daop.tree_offset=(int*) malloc(sizeof(int)*num_trees);
    daop.tree_start=(int*) malloc(sizeof(int)*(num_trees));

     curr_tree=0;
    daop.tree_start[0]=curr_tree;
    daop.tree_offset[0]=0;
    counter=0;
    offset=0;
    for(k=0;k<nrows[0];k++)
      {
	if(daop.tree_matrix[0][k]!=curr_tree)
	  {
	    counter++;
	    curr_tree=daop.tree_matrix[0][k];
	    daop.tree_start[counter]=k;
	    if(curr_tree==0)
	      {
		offset=k;
	      }
	    daop.tree_offset[counter]=offset;
	  }
      }
    
    if(daop.ncat<0){
      if(all_trees[0]==1){
	predict_basic_all(nboost,predictions);
      }else{
	predict_basic(nboost,predictions);
      }
    }

     if(daop.ncat>0){
      if(all_trees[0]==1){
	predict_basic_all_gini(nboost,predictions);
      }else{
	predict_basic_gini(nboost,predictions);
      }


    }

      

 
free(daop.tree_offset);

 
   free_basic();
 
}


double tnode_basic(int tree_indx,int row_number)
{
  
  // For subject in row_number, find terminal node for subject

  int k,tnode,go,j,i,split_var;
  double ff;
  int counter;
  double pred;
  int offset;

  offset=daop.tree_offset[tree_indx];

  tnode=-99;
  k=daop.tree_start[tree_indx];
  go=1;
  counter=0;
   while(go==1)
    {

    
      if(daop.tree_matrix[4][k]<(-10))
	{
	  // terminal node
	  pred=(daop.tree_matrix[3][k]);
	  break;
	}else{
	split_var=((int) daop.tree_matrix[4][k]);

	if((daop.tree_matrix[5][k]<(-10))&&(daop.tree_matrix[14][k]<(-10))) // not a summary variable
	  {
	    if(daop.x[split_var][row_number]<daop.tree_matrix[6][k])
	      {
		// Go Left
		k=((int) daop.tree_matrix[8][k])+offset;
	      }else{
	      // Go Right
	      k=((int) daop.tree_matrix[9][k])+offset;
	    }

	  }else{

	  ff=get_summary_basic2(row_number,(int) daop.tree_matrix[13][k],(int) daop.tree_matrix[14][k],(int) daop.tree_matrix[15][k],daop.tree_matrix[5][k]);
	  if(ff<daop.tree_matrix[6][k])
	    {
		// Go Left
		k=((int) daop.tree_matrix[8][k]+offset);
	      }else{
	      // Go Right
	      k=((int) daop.tree_matrix[9][k]+offset);
	    }
	}

      }
    }

   return(pred);
}




void tnode_basic_gini(int tree_indx,int row_number,double *pred)
{
  // For subject in row_number, find terminal node for subject

  int k,tnode,go,j,i,split_var;
  double ff;
  int counter;
  int offset;
  int NCOL_TREE;

  NCOL_TREE=16;

  offset=daop.tree_offset[tree_indx];

  tnode=-99;
  k=daop.tree_start[tree_indx];
  go=1;
  counter=0;
   while(go==1)
    {

    
      if(daop.tree_matrix[4][k]<(-10))
	{
	  for(j=0;j<daop.ncat;j++)
	    pred[j]=(daop.tree_matrix[NCOL_TREE+j][k]);
	  break;
	}else{
	split_var=((int) daop.tree_matrix[4][k]);

	if((daop.tree_matrix[5][k]<(-10))&&(daop.tree_matrix[14][k]<(-10))) // not a summary variable
	  {
	    if(daop.x[split_var][row_number]<daop.tree_matrix[6][k])
	      {
		// Go Left
		k=((int) daop.tree_matrix[8][k])+offset;
	      }else{
	      // Go Right
	      k=((int) daop.tree_matrix[9][k])+offset;
	    }

	  }else{

	  ff=get_summary_basic2(row_number,(int) daop.tree_matrix[13][k],(int) daop.tree_matrix[14][k],(int) daop.tree_matrix[15][k],daop.tree_matrix[5][k]);
	  if(ff<daop.tree_matrix[6][k])// was 5.. but changed 19/5/2018
	    {
		// Go Left
		k=((int) daop.tree_matrix[8][k]+offset);
	      }else{
	      // Go Right
	      k=((int) daop.tree_matrix[9][k]+offset);
	    }
	}

      }
    }

  
}



void predict_basic(int *ntrees,double *pred)
{

  int i,k;

  for(i=0;i<daop.n;i++)
    pred[i]=0;


  // find terminal node
  for(k=0;k<ntrees[0];k++)
    {
      for(i=0;i<daop.n;i++)
	{
	  pred[i]+=tnode_basic(k,i);
	}

    }

}



void predict_basic_all(int *ntrees,double *pred)
{

  // returns predictions from each tree 
  int i,k;

  for(i=0;i<daop.n;i++)
    pred[i]=0;


  // find terminal node
  for(k=0;k<ntrees[0];k++)
    {
      for(i=0;i<daop.n;i++)
	{
	  pred[i+(daop.n*k)]=tnode_basic(k,i);
	}

    }

}


void predict_basic_all_gini(int *ntrees,double *pred)
{

  // returns predictions from each tree 
  int i,k;
  int ncat,n;

  ncat=daop.ncat;
  n=daop.n;
  
  for(i=0;i<(daop.n*daop.ncat*ntrees[0]);i++)
    pred[i]=0;

 
  // find terminal node
  for(k=0;k<ntrees[0];k++)
    {
      for(i=0;i<daop.n;i++)
	{
	  tnode_basic_gini(k,i,&(pred[i*ncat+(n*ncat*k)]));
 
	}

    }

}

void predict_basic_gini(int *ntrees,double *pred)
{

  // returns predictions from each tree 
  int i,k,m;
  int ncat,n;
  double paux[MAXCAT];

  ncat=daop.ncat;
  n=daop.n;
  
  for(i=0;i<(daop.n*daop.ncat);i++)
    pred[i]=0;


  // find terminal node
  for(k=0;k<ntrees[0];k++)
    {
      //  Rprintf("k=%d \n",k);
  for(i=0;i<daop.n;i++)
    {
      tnode_basic_gini(k,i,&(paux[0]));
      for(m=0;m<daop.ncat;m++){
	pred[i*ncat+m]+=paux[m];
      }
    }

    }

}




// .................................................................
// .. Variable importance 



void varimp_predictor(int j,double *pred)
{

  // Randomly permute predictor 'j' and compute oob predictions for 

  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  
  double *x_orig;
  int *perm_indx;
  double *oobv;
  int i,k,m;


  
  
  
  x_orig=daop.daux1;
  perm_indx=daop.iaux1;
  oobv=daop.daux2;

  for(i=0;i<daop.n;i++)
    perm_indx[i]=i;
  
    permute(perm_indx,daop.n);  
 
    //Rprintf(" vi: 1 \n");
  for(i=0;i<daop.n;i++){
    pred[i]=0;
    oobv[i]=0;
    x_orig[i]=daop.x[j][i];
  }

 
  // permute predictor variable ..
  for(i=0;i<daop.n;i++)
    daop.x[j][i]=x_orig[perm_indx[i]];


  regenerate_history_summary();
 
  // loop through trees and form oob predictions with variable permuted 

  for(k=0;k<daop.nboost;k++)
    {
   
      HASH_FIND_INT(oob,&k,s);
      HASH_ITER(hh,s->rows,r,rtemp){
	m=r->row_number;
	pred[m]=pred[m]+tnode_basic(k,m);
	oobv[m]=oobv[m]+1;
     }

    }
 
  for(m=0;m<daop.n;m++)
    {
      if(oobv[m]>.1&&daop.rf==1)
	pred[m]=pred[m]/oobv[m];
    }

  
  // restore permuted predictor
    for(i=0;i<daop.n;i++) 
      daop.x[j][i]=x_orig[i];
    
}


void varimp(int *nperm,double *res)
{
  int j,k,i;
  double *pred;
  double a;

   GetRNGstate();
  pred=(double*) malloc(sizeof(double)*daop.n);


  for(i=0;i<(daop.n*daop.p);i++)
    res[i]=0;
  for(j=0;j<daop.p;j++)
    {
      for(k=0;k<nperm[0];k++){
	varimp_predictor(j,pred);
	for(i=0;i<daop.n;i++)
	    res[daop.n*j+i]+=pred[i];
	  
      }
    }

  a=(double) nperm[0];
  for(i=0;i<(daop.n*daop.p);i++)
    res[i]=res[i]/a;

  free(pred);
   PutRNGstate();
}






void varimp_predictor_gini(int j,double *pred)
{

  // Randomly permute predictor 'j' and compute oob predictions for 

  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  
  double *x_orig;
  int *perm_indx;
  double *oobv;
  int i,k,m,l;
  double *paux;
  
  paux=(double*) malloc(sizeof(double)*daop.ncat);
  
  
  x_orig=daop.daux1;
  perm_indx=daop.iaux1;
  oobv=daop.daux2;

  for(i=0;i<daop.n;i++)
    perm_indx[i]=i;
  
    permute(perm_indx,daop.n);  
 
  for(i=0;i<daop.n;i++){
    pred[i]=0;
    oobv[i]=0;
    x_orig[i]=daop.x[j][i];
  }

  // permute predictor variable ..
  for(i=0;i<daop.n;i++)
    daop.x[j][i]=x_orig[perm_indx[i]];


  regenerate_history_summary();

  // loop through trees and form oob predictions with variable permuted 

  for(k=0;k<daop.nboost;k++)
    {

      HASH_FIND_INT(oob,&k,s);
      HASH_ITER(hh,s->rows,r,rtemp){
	m=r->row_number;
	tnode_basic_gini(k,m,&(paux[0]));
	  for(l=0;l<daop.ncat;l++)
	    pred[m*daop.ncat+l]=pred[m*daop.ncat+l]+paux[l]; //tnode_basic(k,m);
	oobv[m]=oobv[m]+1;
      }

    }

  for(m=0;m<daop.n;m++)
    {
      if(oobv[m]>.1&&daop.rf==1){
	for(l=0;l<daop.ncat;l++)
	pred[m*daop.ncat+l]=pred[m*daop.ncat+l]/oobv[m];

      }
    }

  
  
  // restore permuted predictor
    for(i=0;i<daop.n;i++) 
      daop.x[j][i]=x_orig[i];


    free(paux);
}


void varimp_gini(int *nperm,double *res)
{
  int j,k,i;
  double *pred;
  double a;

   GetRNGstate();
  pred=(double*) malloc(sizeof(double)*daop.n*daop.ncat);

  for(i=0;i<(daop.n*daop.p*daop.ncat);i++)
    res[i]=0;
  for(j=0;j<daop.p;j++)
    {
      for(k=0;k<nperm[0];k++){
	varimp_predictor_gini(j,pred);
	for(i=0;i<(daop.n*daop.ncat);i++)
	    res[daop.n*daop.ncat*j+i]+=pred[i];
	  
      }
    }

  a=(double) nperm[0];
  for(i=0;i<(daop.n*daop.p*daop.ncat);i++)
    res[i]=res[i]/a;

  free(pred);
   PutRNGstate();
}







void ctree_varimp(double *x,int *id,double *time,int *n,int *p,int *yindx,int *nboost,
	   int *concurrent,int *nconcurrent,  
	   int *historic,int *nhistoric,
		  double *delta,double *delta0,int *ndelta,int *type,int *ncat,double *quantiles,int *nquantiles,
		  double *trees,int *nrows,int *oob_vec,int *nperm,int *rf,double *predictions
	   )
{
  int NCOL_TREE;
  int mtry,nsplit,nodesize,dtry;
  int k,i,curr_tree;
  int num_trees,counter,offset;

  
  
  NCOL_TREE=16;
  if(ncat[0]>0)
    NCOL_TREE=NCOL_TREE+ncat[0];
  
  mtry=1;
  nsplit=1;
  nodesize=1;
  dtry=1;
  daop.rf=rf[0];
  daop.qtry=1;
  if(type[0]==2||type[0]==5)
    daop.quantile_sample=(double*) malloc(sizeof(double)*nquantiles[0]);


 
  read_basic(x,id,time,n,p,yindx,&nsplit,nboost,&mtry,&nodesize,concurrent,nconcurrent,historic,nhistoric);
  free(daop.tree_start); // reallocated below
  create_history(delta,delta0,ndelta,type,quantiles,nquantiles,&dtry);
  // read in oob 
  set_oob(oob_vec);

  daop.ncat=ncat[0];

    
  // Copy training data (daop.x[] allocated in 'read_basic')
  for(k=0;k<p[0];k++)
    {
      daop.x[k]=(double*) malloc(sizeof(double)*daop.n);
      for(i=0;i<daop.n;i++)
	daop.x[k][i]=x[(daop.n)*k+i];
    }
  
  // .. trees
  daop.tree_matrix=(double**) malloc(sizeof(double*)*NCOL_TREE);
  for(k=0;k<NCOL_TREE;k++)
    daop.tree_matrix[k]=&(trees[k*nrows[0]]);

 

  curr_tree=0;

  // count number of trees
  num_trees=1;
  for(k=0;k<nrows[0];k++)
    {
      if(daop.tree_matrix[0][k]!=curr_tree)
	{
	  curr_tree=daop.tree_matrix[0][k];
	  num_trees++;
	}
    }

 
  daop.nboost=num_trees;
  daop.tree_offset=(int*) malloc(sizeof(int)*num_trees);
  daop.tree_start=(int*) malloc(sizeof(int)*(num_trees));
 

  curr_tree=0;
  daop.tree_start[0]=curr_tree;
  daop.tree_offset[0]=0;
  counter=0;
  offset=0;
  for(k=0;k<nrows[0];k++)
    {
      if(daop.tree_matrix[0][k]!=curr_tree)
	{
	  counter++;
	  curr_tree=daop.tree_matrix[0][k];
	  daop.tree_start[counter]=k;
	  if(curr_tree==0)
	    {
	      offset=k;
	    }
	  daop.tree_offset[counter]=offset;
	}
    }
    
  if(daop.ncat<0){
    varimp(nperm,predictions);
  }else{
    varimp_gini(nperm,predictions);
  }
	   
    
  for(k=0;k<p[0];k++)
    free(daop.x[k]);
	
 
  free(daop.tree_offset);
  free_basic();
 
}





void update_terminal_nodes_wrapper(int treeindx)
{

 
  if(daop.family==2)
    {
      update_terminal_nodes(treeindx);
    }

  
}




void update_terminal_nodes(int tree_indx)
{
  int k,i,start,go,ni,ti,end;
  int ntnodes;
  double h;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  double *node_predictions;
  int *tnodes;

  tnodes=(int*) malloc(sizeof(int)*daop.n);
  node_predictions=(double*) malloc(sizeof(double)*daop.n);

  
  ntnodes=0;
  start=daop.tree_start[tree_indx];
  end=daop.tree_end[tree_indx];
  
  for(k=start;k<=end;k++)
    {

      ti=(int) daop.tree_matrix[0][k];
      if(ti==tree_indx)
	{

	  if(daop.tree_matrix[4][k]<(-1))
	    {
	      // THIS is a terminal node ....
	      ni=((int) daop.tree_matrix[2][k]);
	      h=node_prediction_wrapper(ni);
	      daop.tree_matrix[3][k]=(daop.lambda*h);
	    }
	  
	}else{
	break;
      }

    }


  free(tnodes);
  free(node_predictions);

 
}

  
double node_prediction_wrapper(int nodeid)
{
  // Called by: update_terminal_nodes(_wrapper)
  double prediction;
  prediction=0;
  if(daop.family==1){
      prediction=node_prediction_uthash(nodeid);
  }

  if(daop.family==2)
    {
      // logistic regression 
      prediction=node_prediction_logistic(nodeid);
    }

  return(prediction);
}



double node_prediction_logistic(int nodeid)
{
 
  int i;
  double a,b;
  double n,pr;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  a=0;
  b=0;
  n=0;

  HASH_FIND_INT(nodes,&nodeid,s);
  if(s==NULL){}else{
    HASH_ITER(hh,s->rows,r,rtemp)
      {
	i=r->row_number;
	if(daop.train[i]==1)
	  {
	    a+=daop.response[i];  // (y-pr)  
	    pr=1/(1+exp(-daop.predictions[i]));
	    b+=(pr*(1-pr));
	    n++;
	  }
       
      }
  }

  if(n<.5)
    n=1;

  
   a=a/(b+.0001);
  return(a);
}


