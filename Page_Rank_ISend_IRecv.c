#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>

unsigned long int no_of_vertices,no_of_edges;
unsigned int MAXCHAR=150;
char *filename="large_graph.dat";
int *mapping_table;

struct Send_Node
{
    int node;
    double weight; 
};

typedef struct Send_Node SN;

struct LL_Node
{
    int g_node;
    int process_rank;
    struct LL_Node *next;
};
typedef struct LL_Node LL;

LL * Insert_end_node(LL *end_point,int Node,int r)
{
    /*if(Node==1)
    {
        printf("HIIIIII\n");
    }*/
    LL *New;
    New=(LL *)calloc(sizeof(LL),1);
    New->g_node=Node;
    New->process_rank=r;
    New->next=NULL;
    if (end_point==NULL)
        return New;
    end_point->next=New;
    return New;
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

void Initalize_Local_Page_Rank_And_OutNode_Outdegree(double *local_page_rank,int *Outdegree,LL **OutNodes,LL **LastNodes_Out,LL **InNodes,LL **LastNodes_In,int start,int end,int rank)
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
    int count=0;
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
        
        
       
        if(source >= start && source<= end)
        {        
            LastNodes_Out[source-start] = Insert_end_node(LastNodes_Out[source-start],destination,mapping_table[destination-1]);
            Outdegree[source-start]++;
            if(OutNodes[source-start]==NULL)
                OutNodes[source-start] = LastNodes_Out[source-start];
            /*if(rank==0 && destination==1 && source==7)
            {
                printf("Source:%d\n",source);
                LL *curr=OutNodes[source-start];
                while(curr!=NULL)
                {
                    if(curr->g_node==destination)
                    {
                        printf("Process Rank:%d\n",curr->process_rank);
                        break;
                    }
                    curr=curr->next;
                }
            } */   
        }  

        if(destination >= start && destination<= end)
        {    
            LastNodes_In[destination-start] = Insert_end_node(LastNodes_In[destination-start],source,mapping_table[source-1]);
            if(InNodes[destination-start]==NULL)
                InNodes[destination-start] = LastNodes_In[destination-start];
        }     
    }
}

void Update_Partial_Page_Rank(int send_counts,SN *send_buffer,int start,double *partial_page_rank)
{
    SN Temp;
    int i;
    for(i=0;i<send_counts;i++)
    {
        Temp=send_buffer[i];
        partial_page_rank[Temp.node-start]+= Temp.weight;       
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
    
    int *my_size,*displ,*local_a,i,j,k;
    my_size = (int *)malloc(nprocs*sizeof(int *));
    displ = (int *)malloc(nprocs*sizeof(int *));
    mapping_table=(int *)calloc(sizeof(int),no_of_vertices);

    Create_Size_Array(my_size,nprocs);    
    Create_Displacement_Array(displ,my_size,nprocs);

    local_a=(int*)calloc(sizeof(int),my_size[rank]);

    for(i=0;i<my_size[rank];i++)
        local_a[i]=rank;

    MPI_Allgatherv(local_a,my_size[rank],MPI_INT,mapping_table,my_size,displ,MPI_INT,MPI_COMM_WORLD);


    int start=displ[rank]+1,end;

    if (rank==nprocs-1)
        end=no_of_vertices;
    else
        end=displ[rank+1];

    int *Outdegree;
    LL **OutNodes,**LastNodes_Out,**InNodes,**LastNodes_In;
    
    double *local_page_rank;

    Outdegree=(int *)calloc(sizeof(int),my_size[rank]);

    local_page_rank=(double *)calloc(sizeof(double),my_size[rank]);


    OutNodes=(LL **)calloc(sizeof(LL *),my_size[rank]);
    LastNodes_Out=(LL **)calloc(sizeof(LL *),my_size[rank]);

    InNodes=(LL **)calloc(sizeof(LL *),my_size[rank]);
    LastNodes_In=(LL **)calloc(sizeof(LL *),my_size[rank]);

    Initalize_Local_Page_Rank_And_OutNode_Outdegree(local_page_rank,Outdegree,OutNodes,LastNodes_Out,InNodes,LastNodes_In,start,end,rank);
    
    int *send_counts,*recv_counts;
    send_counts=(int *)calloc(sizeof(int),nprocs);
    recv_counts=(int *)calloc(sizeof(int),nprocs);

    LL *curr;

    for(i=0;i<my_size[rank];i++)
    {
        curr=OutNodes[i];

        while(curr!=NULL)
        {
            send_counts[curr->process_rank]+=1;
            curr=curr->next;
        }
    }

    for(i=0;i<my_size[rank];i++)
    {
        curr=InNodes[i];

        while(curr!=NULL)
        {
            recv_counts[curr->process_rank]+=1;
            curr=curr->next;
        }
    }

    /*
    for(j=0;j<nprocs;j++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank==j)
        {
            printf("\nProcess : %d\n",j);
            
            for(i=0;i<nprocs;i++)
                printf("Process:%d Send:%d Recev:%d\n",i,send_counts[i],recv_counts[i]);
        }  
    }*/

    MPI_Request *s_req,*r_req;
    MPI_Status *status;
    SN **send_buffer,**recv_buffer;
    int *last_index_send;
    double *partial_page_rank;

    send_buffer=(SN **)calloc(sizeof(SN *),nprocs);
    recv_buffer=(SN **)calloc(sizeof(SN *),nprocs);
    s_req=(MPI_Request *)calloc(sizeof(MPI_Request),nprocs);
    r_req=(MPI_Request *)calloc(sizeof(MPI_Request),nprocs);
    status=(MPI_Status *)calloc(sizeof(MPI_Status),nprocs);
    

    
    for (i=0;i<nprocs;i++)
        send_buffer[i]=(SN *)calloc(sizeof(SN),send_counts[i]);

    for (i=0;i<nprocs;i++)
        recv_buffer[i]=(SN *)calloc(sizeof(SN),recv_counts[i]);

    MPI_Barrier(MPI_COMM_WORLD);

    double weight,local_error,global_error,d=0.85;
    int process_number,count;
    int *flag_array,flag,flag_test;
    
    
    MPI_Datatype SEND;
    int lengths[2]={1,1};
    MPI_Aint displacement[2];
    displacement[0]=offsetof(SN,node);
    displacement[1]=offsetof(SN,weight);
    MPI_Datatype types[2]={MPI_INT,MPI_DOUBLE};
    MPI_Type_create_struct(2,lengths,displacement,types,&SEND);
    MPI_Type_commit(&SEND);
    
    count=0;
    start_time=MPI_Wtime();

    while(1)
    {    
        //MPI_Barrier(MPI_COMM_WORLD);
        local_error=0;
        global_error=0;

        last_index_send=(int *)calloc(sizeof(int),nprocs);

        for(i=0;i<my_size[rank];i++)
        {
            if(Outdegree[i]!=0)
            {
                weight=local_page_rank[i]/Outdegree[i];

                curr=OutNodes[i];
                while(curr!=NULL)
                {
                    
                    SN New;
                    New.node=curr->g_node;
                    New.weight=weight;
                    
                    process_number=curr->process_rank;
                    send_buffer[process_number][last_index_send[process_number]]=New;
                    last_index_send[process_number]++;
                    curr=curr->next;
                }
            }
        }
    
        for(i=0;i<nprocs;i++)
        {
            if(i!=rank)
                MPI_Isend(send_buffer[i],send_counts[i],SEND,i,0,MPI_COMM_WORLD,&s_req[i]);
        }
        
        for(i=0;i<nprocs;i++)
        {
            if(i!=rank)
                MPI_Irecv(recv_buffer[i],recv_counts[i],SEND,i,0,MPI_COMM_WORLD,&r_req[i]);
        }
        //printf("ISend Irecv Over %d rank.\n",rank);
        partial_page_rank=(double *)calloc(sizeof(double),my_size[rank]);

        
        Update_Partial_Page_Rank(send_counts[rank],send_buffer[rank],start,partial_page_rank);

        flag_array=(int *)calloc(sizeof(int),nprocs);
        flag_array[rank]=1;
        
        while(1)
        {
            flag_test=0;
            for(i=0;i<nprocs;i++)
            {
                if(flag_array[i]!=1)
                {
                    flag_test=1;
                    MPI_Test(&r_req[i],&flag,&status[i]);
                    if(flag==1)
                    {
                        flag_array[i]=1;
                        MPI_Wait(&r_req[i],&status[i]);
                        Update_Partial_Page_Rank(recv_counts[i],recv_buffer[i],start,partial_page_rank);
                    }  
                }
            }
            if(flag_test==0)
                break;
        }
        
        /*
        for(i=0;i<nprocs;i++)
        {
            if(rank!=i)
                Update_Partial_Page_Rank(recv_counts[i],recv_buffer[i],start,partial_page_rank);

        }*/
        double ans;

        for(i=0;i<my_size[rank];i++)
        {
            ans=((1-d)/(double)no_of_vertices)+(d*partial_page_rank[i]);
            local_error+=fabs(local_page_rank[i]-ans);
            local_page_rank[i]=ans;
        }

        free(partial_page_rank);
        free(last_index_send);

        MPI_Allreduce(&local_error,&global_error,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        if(global_error < 1e-4)
            break;
        else
        {
            count++;
            /*if(rank==0)
                printf("%d Itration %g Error\n",count,global_error);*/
        }
        if(count > 35)
            break;
        
    }
    
    my_time=MPI_Wtime()-start_time;
    MPI_Reduce(&my_time,&exec_time,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if(rank==0)
    {
        
        /*for(i=0;i<my_size[rank];i++)
            printf("%d node %g rank\n",i+1,local_page_rank[i]);*/
        printf("Total Time:%g sec",exec_time/nprocs);    
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}