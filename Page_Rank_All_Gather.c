#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

unsigned long int no_of_vertices,no_of_edges;
unsigned int MAXCHAR=150;
char *filename="graph.dat";

struct LL_Node
{
    int g_node;
    struct LL_Node *next;
};
typedef struct LL_Node LL;

LL * Insert_end(LL *end_point,int Node)
{
    LL *New;
    New=(LL *)calloc(sizeof(LL),1);
    New->g_node=Node;
    New->next=NULL;
    if (end_point==NULL)
        return New;
    end_point->next=New;
    return end_point;
}

int isNum(char c)
{
    if (c=='0' || c=='1' ||c=='2' ||c=='3' ||c=='4' ||c=='5' ||c=='6' ||c=='7' ||c=='8' ||c=='9')
        return 1;
    else 
        return 0;
}

char *strrev(char *str)
{
      char *p1, *p2;

      if (! str || ! *str)
            return str;
      for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
      {
            *p1 ^= *p2;
            *p2 ^= *p1;
            *p1 ^= *p2;
      }
      return str;
}

int Read_Last_Int_From_Line(char *str)
{
    char Vertices_String[MAXCHAR];
    unsigned int len,i;
    len=strlen(str);
    Vertices_String[0]='\0';
    i=len-2;

    while (i >= 0)
    {
        char c=str[i];
        if(isNum(c))
        {
            strncat(Vertices_String,&c,1);
            i-=1;
        }
        else
            break;
    }
    strrev(Vertices_String);
    return atoi(Vertices_String);
}

void Initalize_No_Of_Vertices_Edges()
{
    FILE *GF;
        char str[MAXCHAR];
        unsigned int len,i;
        

        GF=fopen(filename,"r");

        fgets(str,MAXCHAR,GF);
        //First Line is useless 
        fgets(str,MAXCHAR,GF);
        //Read Number Of Vertices  
        no_of_vertices=Read_Last_Int_From_Line(str);   
        
        fgets(str,MAXCHAR,GF);
        //Read Number Of Edges  
        no_of_edges=Read_Last_Int_From_Line(str);

        //printf("Number Of Vertices:%lu\n",no_of_vertices);
        //printf("Number Of Edges:%lu\n",no_of_edges);

        fclose(GF);
}

void Create_Size_Array(int *my_size,int nprocs)
{
    int num_element = no_of_vertices/nprocs,i;
        for(i=0; i<nprocs; i++)
            my_size[i] = num_element;

        int rem = no_of_vertices%nprocs;
        i=0;
        while(rem--){
            my_size[i]++;
            i++;
        }
}

void Create_Displacement_Array(int *displ,int *my_size,int nprocs)
{
    int i;
    displ[0] = 0;
    for(i=1; i<nprocs; i++)
        displ[i] = displ[i-1] + my_size[i-1];
}

void Initalize_Local_Page_Rank_And_OutNode_Outdegree(double *local_page_rank,int *Outdegree,LL **OutNodes,LL **LastNodes,int start,int end)
{
    FILE *GF;
    char str[MAXCHAR];
    unsigned int len,i;

    for(i=0;i<end-start+1;i++)
        local_page_rank[i]=1/(double)no_of_vertices;

    GF=fopen(filename,"r");

    char *token;
    int source,destination,weight;
    char *info[4];
    
    while (fgets(str, MAXCHAR, GF) != NULL)
    {
        if(str[0]!='a')
            continue;
        else
        {
            i=0;
            token=strtok(str," ");
            while(token!=NULL)
            {
                info[i]=token;
                i++;
                token = strtok(NULL, " ");
            }
        }
        
        info[3][strlen(info[3])-1]='\0';
        source=atoi(info[1]);
        destination=atoi(info[2]);
        weight=atoi(info[3]);
        
        /*
        local_page_rank[source]=1/(double)no_of_vertices;
        local_page_rank[destination]=1/(double)no_of_vertices;
        */
        if(destination >= start && destination<= end)
        {    
            LastNodes[destination-start] = Insert_end(LastNodes[destination-start],source);
            if(OutNodes[destination-start]==NULL)
                OutNodes[destination-start] = LastNodes[destination-start];
        }
        if(source >= start && source<=end)
            Outdegree[source-start]++;     
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);
    
    int rank,nprocs;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    double start_time,my_time,exec_time;
    

    Initalize_No_Of_Vertices_Edges();
    
    int *my_size,*displ,i,j,k;
    my_size = (int *)malloc(nprocs*sizeof(int *));
    displ = (int *)malloc(nprocs*sizeof(int *));

    Create_Size_Array(my_size,nprocs);    
    Create_Displacement_Array(displ,my_size,nprocs);
        
    int start=displ[rank]+1,end;

    if (rank==nprocs-1)
        end=no_of_vertices;
    else
        end=displ[rank+1];

    int *Outdegree;
    LL **OutNodes,**LastNodes;
    
    double *local_page_rank;

    Outdegree=(int *)calloc(sizeof(int),my_size[rank]);

    local_page_rank=(double *)calloc(sizeof(double),no_of_vertices);


    OutNodes=(LL **)calloc(sizeof(LL *),my_size[rank]);
    LastNodes=(LL **)calloc(sizeof(LL *),my_size[rank]);


    Initalize_Local_Page_Rank_And_OutNode_Outdegree(local_page_rank,Outdegree,OutNodes,LastNodes,start,end);
    
    MPI_Barrier(MPI_COMM_WORLD);

    double *Own_Part,*global_weight,sum,global_sum;
    double d=0.85,local_error,global_error;
    int count=0,index;

    start_time=MPI_Wtime();
    Own_Part=(double *)calloc(sizeof(double),my_size[i]);
    global_weight=(double *)calloc(sizeof(double),no_of_vertices);  
   
    LL *curr;
    
    while(1)
    {  

        local_error=0;
        global_error=0;
  
        for(i=0;i<my_size[rank];i++)
        {
            if(Outdegree[i]!=0)
                Own_Part[i]=local_page_rank[i]/Outdegree[i];
            else
                Own_Part[i]=0;
        }
        MPI_Allgatherv(Own_Part,my_size[rank],MPI_DOUBLE,global_weight,my_size,displ,MPI_DOUBLE,MPI_COMM_WORLD);

        for(i=0;i<my_size[rank];i++)
        {
            sum=0;
            curr=OutNodes[i];

            while (curr!=NULL)
            {
                index=curr->g_node;
                sum+=global_weight[index-1];
                curr=curr->next;
            }

            sum=((1-d)/no_of_vertices)+ d*sum;
            local_error+=fabs(sum-local_page_rank[i]);
            local_page_rank[i]=sum;
        }

        MPI_Allreduce(&local_error,&global_error,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        if(global_error < 1e-10)
            break;
        count++;
        if(count<35)
        {
            if (rank==0)
                printf("%d Itration %g error\n",count,global_error);
        }
        else
            break;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    my_time=MPI_Wtime()-start_time;
    MPI_Reduce(&my_time, &exec_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if(rank==0)
    {
        /*
        for(i=0;i<my_size[rank];i++)
        {
            printf("%g\n",local_page_rank[i]);
        }
        */
        //printf("Sum:%g\n",global_sum);        
        printf("Time:%g sec\n",exec_time/nprocs);    
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;

}