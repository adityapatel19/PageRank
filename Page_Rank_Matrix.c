#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

unsigned long int no_of_vertices,no_of_edges;
unsigned int MAXCHAR=150;
char *filename="graph.dat";
int *displ; 

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

int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);
    
    int rank,nprocs;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    displ=(int *)malloc(sizeof(int)*(nprocs+1));
    double start_time,my_time,exec_time;
    start_time=MPI_Wtime();
   
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

        //fclose(GF);
        int *my_size;
        my_size = (int *)malloc(nprocs*sizeof(int *));

        int num_element = no_of_vertices/nprocs;
        for(int i=0; i<nprocs; i++)
            my_size[i] = num_element;

        int rem = no_of_vertices%nprocs;
        i=0;
        while(rem--){
            my_size[i]++;
            i++;
        }

        displ[0] = 0;
        for(i=1; i<nprocs; i++)
            displ[i] = displ[i-1] + my_size[i-1];
        displ[nprocs]=no_of_vertices;    

    //MPI_Barrier(MPI_COMM_WORLD);
    
    int start=displ[rank]+1,end=displ[rank+1];
    
    double **weight_matrix;

    weight_matrix=(double **)calloc(sizeof(double *),(no_of_vertices));

    for (int i=0;i<no_of_vertices;i++)
        weight_matrix[i]=(double *)calloc(sizeof(double),(my_size[rank]));

    
    //GF=fopen(filename,"r");
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
        
        if(source >= start && source<= end)
        {
           weight_matrix[destination-1][source-start]=weight;   
        }     
    }   
    
   
    for (int i=0;i<my_size[rank];i++)
    {
        unsigned long int sum=0;
        for (int j=0;j<no_of_vertices;j++)
            sum+=weight_matrix[j][i];
        if (sum!=0)
        {
            for(int j=0;j<no_of_vertices;j++)
                weight_matrix[j][i]=weight_matrix[j][i] / sum;
        }
    }
    
    double *my_p,*new_p,*local_p,*global_p,local_diff,global_diff;
    my_p=(double *)calloc(sizeof(double),my_size[rank]);
    local_p=(double *)calloc(sizeof(double),no_of_vertices);
    global_p=(double *)calloc(sizeof(double),no_of_vertices);
    new_p=(double *)calloc(sizeof(double),my_size[rank]);

    

    

    double d=1/(double)nprocs;

    for(int i=0;i<my_size[rank];i++)
        my_p[i]= d;
    
    unsigned long int count=0;
    while(1)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        
        for(i=0;i<no_of_vertices;i++)
        {
            local_p[i]=0;
            for(int j=0;j<my_size[rank];j++)
                if (weight_matrix[i][j]!=0)
                    local_p[i]+=weight_matrix[i][j]*my_p[j];    
        }
        
        MPI_Reduce(local_p,global_p,no_of_vertices,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        
        MPI_Scatterv(global_p,my_size,displ,MPI_DOUBLE,new_p,my_size[rank],MPI_DOUBLE,0,MPI_COMM_WORLD);
        
        /*
        local_diff=0;
        for(i=0;i<my_size[rank];i++)
            local_diff+=fabs(my_p[i]-new_p[i]);
        */
        MPI_Allreduce(&local_diff,&global_diff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        /*
        if(global_diff < 1e-1)
            break;
        else
        {
            count++;
            for(i=0;i<my_size[rank];i++)
                my_p[i]=new_p[i];
            if(count > 5)
                break;
        }*/
        count++;
            for(i=0;i<my_size[rank];i++)
                my_p[i]=new_p[i];
            if(count > 0)
                break;
    }
    my_time=MPI_Wtime()-start_time;
    double local_sum=0,global_sum;;

    for(i=0;i<my_size[rank];i++)
        local_sum+=new_p[i];

    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Reduce(&local_sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    //my_time=MPI_Wtime()-start_time;
    MPI_Reduce(&my_time, &exec_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank==0)
    {
        /*
        for(i=0;i<no_of_vertices;i++)
            printf("%g\n",global_p[i]);*/
        
        printf("Total Time=%g msec\n",1000*(exec_time/nprocs));
        //printf("Sum:%g\n",global_sum);
        //printf("Total No Of Itration:%lu\n",count);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;

}