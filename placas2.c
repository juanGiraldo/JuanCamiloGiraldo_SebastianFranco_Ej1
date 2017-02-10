#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* Se definen las constantes del problema */
#define L 5.0
#define l 2.0
#define d 1.0
#define h 0.02
#define V_o 100.0

int matrix(int fila, int columna, int n_dist);

int main(int argc, char** argv){

  MPI_Init(NULL, NULL);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  /* Se inicializan las variables y arreglos que se van a usar en el problema */
  /* N es el numero de iteraciones del metodo de relajacion */
  int N = (int)(2*L*L/(h*h));
  /* n_dis representa el numero de filas y de columnas de la matriz de discretizacion */
  int n_dis = (int)((L/h)+1);

  /* V es el arreglo que representa al potencial */
  float *V;
  V = malloc((n_dis*n_dis)*sizeof(float));
  /* V_temp es un arreglo temporal para cambiar entre iteraciones */
  float *V_temp;
  V_temp= malloc((n_dis*n_dis)*sizeof(float));
  /* verif es un arreglo que funciona como un arreglo de booleanos, con el fin de decidir cuales casillas se evaluan y cuales no */
  int *verif ;
  verif= malloc((n_dis*n_dis)*sizeof(int));

  float *front_down;
  front_down = malloc(n_dis*sizeof(float));

  /* Ciclos que llenan la matriz inicial */
  /* i es una variable auxiliar para los ciclos */
  int i;
  int j;
  /* 'for' que llena los arreglos con ceros */
  for(i=0; i<(n_dis*n_dis); i++){
    V[i] = 0.0;
    V_temp[i] = 0.0;
    verif[i] = 0;
  }
  /* 'for' que inicializa con 1 los bordes de la matriz */
  for(i=0; i<n_dis; i++){
    verif[matrix(0, i, n_dis)] = 1;
    verif[matrix(n_dis-1, i, n_dis)] = 1;
    verif[matrix(i, 0, n_dis)] = 1;
    verif[matrix(i, n_dis-1, n_dis)] = 1;
  }
  /* 'for' que inicializa con V/2 las placas */
  for(i=(int)(L-l)/(2*h); i<=(int)(((L-l)/(2*h))+(l/h)); i++){
    V[matrix((int)((L-d)/(2.0*h)), i, n_dis)] = V_o/2.0;
    V[matrix((int)(((L-d)/(2.0*h))+(d/h)), i, n_dis)] = -V_o/2.0;
    V_temp[matrix((int)((L-d)/(2.0*h)), i, n_dis)] = V_o/2.0;
    V_temp[matrix((int)(((L-d)/(2.0*h))+(d/h)), i, n_dis)] = -V_o/2.0;
    verif[matrix((int)((L-d)/(2.0*h)), i, n_dis)] = 1;
    verif[matrix((int)(((L-d)/(2.0*h))+(d/h)), i, n_dis)] = 1;
  }

  float *Vrank;
  Vrank=malloc(n_dis*(n_dis/4+1)*sizeof(float));
  float *Vrank_temp;
  Vrank_temp=malloc(n_dis*(n_dis/4+1)*sizeof(float));

  if(world_rank<world_size-1){
    for(i=0;i<n_dis/4+1;i++){
      for(j=0;j<n_dis;j++){
        Vrank[matrix(i,j,n_dis)]=V[matrix(world_rank*n_dis/4+i,j,n_dis)];
      }
    }
  }
  else if(world_rank==world_size-1){
    for(i=0;i<n_dis/4;i++){
      for(j=0;j<n_dis;j++){
        Vrank[matrix(i,j,n_dis)]=V[matrix(world_rank*n_dis/4+i,j,n_dis)];
      }
    }
  }

  /* Ciclos que calculan el potencial */
  /* j y k son variables auxiliares para los ciclos */
  int k;
  for(i=0; i<N; i++){
    if(world_rank!=world_size-1){
      for(j=1; j<n_dis/4; j++){
        for(k=1;k<n_dis-1; k++){
        	if(Vrank[matrix(j, k, n_dis)]!=V_o/2 && Vrank[matrix(j, k, n_dis)]!=-V_o/2){
        	  Vrank_temp[matrix(j, k, n_dis)] = (Vrank[matrix(j, k-1, n_dis)]+Vrank[matrix(j, k+1, n_dis)]+Vrank[matrix(j-1, k, n_dis)]+Vrank[matrix(j+1, k, n_dis)])/4;
        	}
        }
      }
      if (world_rank==0){
        MPI_Recv(&front_down, 1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(j=0;j<n_dis;j++){
            MPI_Recv(&front_down[j], 1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for(j=1; j<n_dis-1; j++){
          Vrank_temp[matrix(n_dis/4,j,n_dis)]=Vrank[matrix(n_dis/4-1,j,n_dis)]+Vrank[matrix(n_dis/4,j-1, n_dis)]+Vrank[matrix(n_dis/4,j+1, n_dis)]+front_down[j];
        }
        for(j=0;j<n_dis;j++){
          MPI_Send(&Vrank_temp, 1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
        }
      }
      if (world_rank==1){
        for(j=0;j<n_dis;j++){
          MPI_Send(&Vrank_temp[matrix(1,j,n_dis)], 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        }
        for(j=0;j<n_dis;j++){
          MPI_Recv(&front_down[j], 1, MPI_FLOAT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for(j=0;j<n_dis;j++){
          Vrank_temp[matrix(n_dis/4, j, n_dis)]=Vrank[matrix(n_dis/4-1,j,n_dis)]+Vrank[matrix(n_dis/4,j-1, n_dis)]+Vrank[matrix(n_dis/4,j+1, n_dis)]+front_down[j];
        }
        for(j=0;j<n_dis;j++){
          MPI_Recv(&Vrank_temp[matrix(0,j,n_dis)], 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for(j=0;j<n_dis;j++){
          MPI_Send(&Vrank_temp[matrix(n_dis/4,j,n_dis)], 1, MPI_FLOAT, 2, 0, MPI_COMM_WORLD);
        }
      }
      if(world_rank==2){
        for(j=0;j<n_dis;j++){
          MPI_Send(&Vrank_temp[matrix(1,j,n_dis)], 1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
        }
        for(j=0;j<n_dis;j++){
          MPI_Recv(&front_down[j], 1, MPI_FLOAT, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for(j=0;j<n_dis;j++){
          Vrank_temp[matrix(n_dis/4, j, n_dis)]=Vrank[matrix(n_dis/4-1,j,n_dis)]+Vrank[matrix(n_dis/4,j-1, n_dis)]+Vrank[matrix(n_dis/4,j+1, n_dis)]+front_down[j];
        }
        for(j=0;j<n_dis;j++){
          MPI_Recv(&Vrank_temp[matrix(0,j,n_dis)], 1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        for(j=0;j<n_dis;j++){
          MPI_Send(&Vrank_temp[matrix(n_dis/4,j,n_dis)], 1, MPI_FLOAT, 3, 0, MPI_COMM_WORLD);
        }
      }
    }
    else if(world_rank==world_size-1){
      for(j=1; j<n_dis/4-1; j++){
        for(k=1;n_dis-1; k++){
          if(Vrank[matrix(j, k, n_dis)]!=V_o/2 && Vrank[matrix(j, k, n_dis)]!=-V_o/2){
            Vrank_temp[matrix(j, k, n_dis)] = (Vrank[matrix(j, k-1, n_dis)]+Vrank[matrix(j, k+1, n_dis)]+Vrank[matrix(j-1, k, n_dis)]+Vrank[matrix(j+1, k, n_dis)])/4;
          }
        }
      }
      for(j=0;j<n_dis;j++){
        MPI_Send(&Vrank_temp[matrix(1,j,n_dis)], 1, MPI_FLOAT, 2, 0, MPI_COMM_WORLD);
      }
      for(j=0;j<n_dis;j++){
        MPI_Recv(&Vrank_temp[matrix(0,j,n_dis)], 1, MPI_FLOAT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
    for(j=0; j<(n_dis*(n_dis/4+1)); j++){
      Vrank[j] = Vrank_temp[j];
    }
  }

  if(world_rank!=0){
    for(i=0;i<n_dis/4;i++){
      for(j=0;j<n_dis;j++){
        MPI_Send(&Vrank[matrix(i,j,n_dis)], 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
      }
    }
  }

    if(world_rank==0){
      for(i=0;i<n_dis/4;i++){
        for(j=0;j<n_dis;j++){
          V[matrix(i,j,n_dis)]=Vrank[matrix(i,j,n_dis)];
          MPI_Recv(&V[matrix(n_dis/4+i,j,n_dis)], 1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(&V[matrix(n_dis/2+i,j,n_dis)], 1, MPI_FLOAT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(&V[matrix(3*n_dis/4+i,j,n_dis)], 1, MPI_FLOAT, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      }
    }
if(world_rank==0){
  float dVx, dVy;
  for(i=0; i<n_dis; i++){
    for(j=0; j<n_dis; j++){
      if(i==0){
	dVy = -(V[matrix(i+1, j, n_dis)]-V[matrix(i, j, n_dis)])/h;
      }
      else if(i==(n_dis-1)){
	dVy = -(V[matrix(i, j, n_dis)]-V[matrix(i-1, j, n_dis)])/h;
      }
      else{
	dVy = -(V[matrix(i+1, j, n_dis)]-V[matrix(i-1, j, n_dis)])/(2*h);
      }
      if(j==0){
	dVx = -(V[matrix(i, j+1, n_dis)]-V[matrix(i, j, n_dis)])/h;
      }
      else if(j==(n_dis-1)){
	dVx = -(V[matrix(i, j, n_dis)]-V[matrix(i, j-1, n_dis)])/h;
      }
      else{
	dVx = -(V[matrix(i, j+1, n_dis)]-V[matrix(i, j-1, n_dis)])/(2*h);
      }
      printf("%f %f %f\n", V[matrix(i, j, n_dis)], dVx, dVy);
    }
  }
}

  MPI_Finalize();

  return 0;

}

int matrix(int fila, int columna, int n_dist){
  return (n_dist*fila)+columna;
}
