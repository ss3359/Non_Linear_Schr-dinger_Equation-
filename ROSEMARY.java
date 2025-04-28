import java.util.*;
import java.math.*;
import javax.swing.*;
import java.io.*;
import java.sql.*;


class TimeIndependentSchrodinger{

    double t=0.1,x=0.8,y=0.8,z=0.8;

    double h=6.67e-34; //Plank Constant in Joule seconds. 
    double m=1.67e-27; //Mass of the particle (electrons) kgs.
    int size=100; //Number of Iterations in RK Algorithm. 
    double dt=0.01,dx=0.01,dy=0.01,dz=0.01; //Time step;
    double h_step=0.1;
    double kappa=1e-3; //Spring Constant (N/m)
    double omega=Math.sqrt(kappa/m);

    //Real Component in the PSI function in terms of x, y, and z
    double[][][][] u_x=new double[2][size][size][size]; 
    double[][][][] u_y=new double[2][size][size][size]; 
    double[][][][] u_z=new double[2][size][size][size]; 



  //Imaginary Component of PSI function in terms of x,y,and z
    double[][][][] v_x=new double[2][size][size][size]; 
    double[][][][] v_y=new double[2][size][size][size]; 
    double[][][][] v_z=new double[2][size][size][size]; 

    //Potential Energy Function
    public double V(int i,int j,int k,int l){
        double new_x=x+(i*dx),new_y=y+(j*dy),new_z=z+(k*dz) ;
        double length_squared=(new_x*new_x)+(new_y*new_y)+(new_z*new_z);
        return (0.5)*m*(omega*omega)*(length_squared);
    }

    //First Spatial Derivatives
    public double ddx(double temp[][][][],int i, int j, int k, int l){
        if(j+1>=size|| j-1<0){
            return 0;
        }
        return (temp[i][j+1][k][l]-temp[i][j-1][k][l])/(2*dt);
    }
    public double ddy(double temp[][][][],int i, int j, int k, int l){
        if(k+1>=size|| k-1<0){
            return 0;
        }
        return (temp[i][j][k+1][l]-temp[i][j][k-1][l])/(2*dt);
    }
    public double ddz(double temp[][][][],int i, int j, int k, int l){
        if(l+1>=size|| l-1<0){
            return 0;
        }
        return (temp[i][j][k][l+1]-temp[i][j][k][l-1])/(2*dt);
    }

    //Second Spatial Derivatives
    public double d2dx2(double temp[][][][],int i, int j, int k, int l){
        if(j+1>=size|| j-1<0){
            return 0;
        }
        return (temp[i][j+1][k][l]-2*temp[i][j][k][l]+temp[i][j-1][k][l])/(dt*dt);
    }
    public double d2dy2(double temp[][][][],int i, int j, int k, int l){
        if(k+1>=size|| k-1<0){
            return 0;
        }
        return (temp[i][j][k+1][l]-2*temp[i][j][k][l]+temp[i][j][k-1][l])/(dt*dt);
    }
    public double d2dz2(double temp[][][][],int i, int j, int k, int l){
        if(l+1>=size|| l-1<0){
            return 0;
        }
        return (temp[i][j][k][l+1]-2*temp[i][j][k][l]+temp[i][j][k][l-1])/(dt*dt);
    }

    public double Laplacian(double temp[][][][],int i, int j, int k, int l){
        return d2dx2(temp,i,j,k,l)+d2dy2(temp,i,j,k,l)+d2dz2(temp,i,j,k,l);
    }

    public void FillGridConditions(double temp1[][][][], double[][][][]temp2){
        for(int i=0; i<2; i++){
            for(int j=0; j<size; j++){
                for(int k=0; k<size; k++){
                    for(int l=0; l<size; l++){
                        if(i==0){
                            double new_x=x+(i*dx),new_y=y+(j*dy),new_z=z+(k*dz) ;
                            temp1[i][j][k][l]=Math.exp(-((new_x*new_x)+(new_y*new_y)+(new_z*new_z)));
                            temp2[i][j][k][l]=Math.exp(-((new_x*new_x)+(new_y*new_y)+(new_z*new_z)));    
                        }
                        else{
                            temp1[i][j][k][l]=0.0;
                            temp2[i][j][k][l]=0.0;
                        }
                    }
                }
            }
        }

    }

    //The Method To Update Points On The u and v grids. This will apply the Runge-Kutta Method
    public void UpdateGrids(){
    
        //The Real Component variables
        //Step sizes for x component
        double [][][][] J1_u=new double[2][size][size][size];
        double [][][][] J2_u=new double[2][size][size][size];
        double [][][][] J3_u=new double[2][size][size][size];
        double [][][][] J4_u=new double[2][size][size][size];  
        //Step sizes for y component
        double [][][][] K1_u=new double[2][size][size][size];
        double [][][][] K2_u=new double[2][size][size][size];
        double [][][][] K3_u=new double[2][size][size][size];
        double [][][][] K4_u=new double[2][size][size][size];  

        //Step sizes for z component
        double [][][][] L1_u=new double[2][size][size][size];
        double [][][][] L2_u=new double[2][size][size][size];
        double [][][][] L3_u=new double[2][size][size][size];
        double [][][][] L4_u=new double[2][size][size][size];  

        //The Imaginary Component variables
        //Step sizes for x component
        double [][][][] J1_v=new double[2][size][size][size];
        double [][][][] J2_v=new double[2][size][size][size];
        double [][][][] J3_v=new double[2][size][size][size];
        double [][][][] J4_v=new double[2][size][size][size];  
        //Step sizes for y component
        double [][][][] K1_v=new double[2][size][size][size];
        double [][][][] K2_v=new double[2][size][size][size];
        double [][][][] K3_v=new double[2][size][size][size];
        double [][][][] K4_v=new double[2][size][size][size];  

        //Step sizes for z component
        double [][][][] L1_v=new double[2][size][size][size];
        double [][][][] L2_v=new double[2][size][size][size];
        double [][][][] L3_v=new double[2][size][size][size];
        double [][][][] L4_v=new double[2][size][size][size];  

        //Temp Arrays for x, y, and z 
        double[][][][] temp_x_u=new double[2][size][size][size];
        double[][][][] temp_y_u=new double[2][size][size][size];
        double[][][][] temp_z_u=new double[2][size][size][size];

        double[][][][] temp_x_v=new double[2][size][size][size];
        double[][][][] temp_y_v=new double[2][size][size][size];
        double[][][][] temp_z_v=new double[2][size][size][size];

        //First, we will apply the initial and boundary conditions for the variables x, y, and z in the real 
        //and imaginary components of the PSI function. 
        FillGridConditions(u_x,v_x);
        FillGridConditions(u_y,v_y);
        FillGridConditions(u_z,v_z);

        


        //We will update the points on the grid using the Runge-Kutta Method
        //First Step
        for(int i=0; i<2; i++){
            for(int j=0; j<size; j++){
                for(int k=0; k<size; k++){
                    for(int l=0; l<size; l++){
                        J1_u[i][j][k][l]=(-h/(2*m))*Laplacian(v_x,i,j,k,l)+V(i,j,k,l)*v_x[i][j][k][l];
                        K1_u[i][j][k][l]=(-h/(2*m))*Laplacian(v_y,i,j,k,l)+V(i,j,k,l)*v_y[i][j][k][l];
                        L1_u[i][j][k][l]=(-h/(2*m))*Laplacian(v_z,i,j,k,l)+V(i,j,k,l)*v_z[i][j][k][l];

                        J1_v[i][j][k][l]=(h/(2*m))*Laplacian(u_x,i,j,k,l)+V(i,j,k,l)*u_x[i][j][k][l];
                        K1_v[i][j][k][l]=(h/(2*m))*Laplacian(u_y,i,j,k,l)+V(i,j,k,l)*u_y[i][j][k][l];
                        L1_v[i][j][k][l]=(h/(2*m))*Laplacian(u_z,i,j,k,l)+V(i,j,k,l)*u_z[i][j][k][l];
                    }
                }                            
             }
        }


        //Second Step
        for(int i=0; i<2; i++){
            for(int j=0; j<size; j++){
                for(int k=0; k<size; k++){
                    for(int l=0; l<size; l++){
                       
                        temp_x_u[i][j][k][l]=u_x[i][j][k][l]+(dx/2)*J1_u[i][j][k][l];    
                        temp_y_u[i][j][k][l]=u_y[i][j][k][l]+(dy/2)*K1_u[i][j][k][l];
                        temp_z_u[i][j][k][l]=u_z[i][j][k][l]+(dz/2)*L1_u[i][j][k][l];

                        temp_x_v[i][j][k][l]=v_x[i][j][k][l]+(dx/2)*J1_v[i][j][k][l];    
                        temp_y_v[i][j][k][l]=v_y[i][j][k][l]+(dy/2)*K1_v[i][j][k][l];
                        temp_z_v[i][j][k][l]=v_z[i][j][k][l]+(dz/2)*L1_v[i][j][k][l];

                        J2_u[i][j][k][l]=(-h/(2*m))*Laplacian(temp_x_u,i,j,k,l)+V(i,j,k,l)*temp_x_u[i][j][k][l];
                        K2_u[i][j][k][l]=(-h/(2*m))*Laplacian(temp_y_u,i,j,k,l)+V(i,j,k,l)*temp_y_u[i][j][k][l];
                        L2_u[i][j][k][l]=(-h/(2*m))*Laplacian(temp_z_u,i,j,k,l)+V(i,j,k,l)*temp_z_u[i][j][k][l];

                        J2_v[i][j][k][l]=(h/(2*m))*Laplacian(temp_x_v,i,j,k,l)+V(i,j,k,l)*temp_x_v[i][j][k][l];
                        K2_v[i][j][k][l]=(h/(2*m))*Laplacian(temp_y_v,i,j,k,l)+V(i,j,k,l)*temp_y_v[i][j][k][l];
                        L2_v[i][j][k][l]=(h/(2*m))*Laplacian(temp_z_v,i,j,k,l)+V(i,j,k,l)*temp_z_v[i][j][k][l];
                    }
                }                            
             }
        }
             //Third Step 
             for(int i=0; i<2; i++){
                for(int j=0; j<size; j++){
                    for(int k=0; k<size; k++){
                        for(int l=0; l<size; l++){
                            
    
                            temp_x_u[i][j][k][l]=u_x[i][j][k][l]+(dx/2)*J2_u[i][j][k][l];    
                            temp_y_u[i][j][k][l]=u_y[i][j][k][l]+(dy/2)*K2_u[i][j][k][l];
                            temp_z_u[i][j][k][l]=u_z[i][j][k][l]+(dz/2)*L2_u[i][j][k][l];
    
                            temp_x_v[i][j][k][l]=v_x[i][j][k][l]+(dx/2)*J2_v[i][j][k][l];    
                            temp_y_v[i][j][k][l]=v_y[i][j][k][l]+(dy/2)*K2_v[i][j][k][l];
                            temp_z_v[i][j][k][l]=v_z[i][j][k][l]+(dz/2)*L2_v[i][j][k][l];
    
                            J3_u[i][j][k][l]=(-h/(2*m))*Laplacian(temp_x_u,i,j,k,l)+V(i,j,k,l)*temp_x_u[i][j][k][l];
                            K3_u[i][j][k][l]=(-h/(2*m))*Laplacian(temp_y_u,i,j,k,l)+V(i,j,k,l)*temp_y_u[i][j][k][l];
                            L3_u[i][j][k][l]=(-h/(2*m))*Laplacian(temp_z_u,i,j,k,l)+V(i,j,k,l)*temp_z_u[i][j][k][l];
    
                            J3_v[i][j][k][l]=(h/(2*m))*Laplacian(temp_x_v,i,j,k,l)+V(i,j,k,l)*temp_x_v[i][j][k][l];
                            K3_v[i][j][k][l]=(h/(2*m))*Laplacian(temp_y_v,i,j,k,l)+V(i,j,k,l)*temp_y_v[i][j][k][l];
                            L3_v[i][j][k][l]=(h/(2*m))*Laplacian(temp_z_v,i,j,k,l)+V(i,j,k,l)*temp_z_v[i][j][k][l];
                        }
                    }                            
                 }
                
        }

        //Fourth Step
        for(int i=0; i<2; i++){
            for(int j=0; j<size; j++){
                for(int k=0; k<size; k++){
                    for(int l=0; l<size; l++){
                       
                        temp_x_u[i][j][k][l]=u_x[i][j][k][l]+(dx)*J3_u[i][j][k][l];    
                        temp_y_u[i][j][k][l]=u_y[i][j][k][l]+(dy)*K3_u[i][j][k][l];
                        temp_z_u[i][j][k][l]=u_z[i][j][k][l]+(dz)*L3_u[i][j][k][l];

                        temp_x_v[i][j][k][l]=v_x[i][j][k][l]+(dx)*J3_v[i][j][k][l];    
                        temp_y_v[i][j][k][l]=v_y[i][j][k][l]+(dy)*K3_v[i][j][k][l];
                        temp_z_v[i][j][k][l]=v_z[i][j][k][l]+(dz)*L3_v[i][j][k][l];

                        J4_u[i][j][k][l]=(-h/(2*m))*Laplacian(temp_x_u,i,j,k,l)+V(i,j,k,l)*temp_x_u[i][j][k][l];
                        K4_u[i][j][k][l]=(-h/(2*m))*Laplacian(temp_y_u,i,j,k,l)+V(i,j,k,l)*temp_y_u[i][j][k][l];
                        L4_u[i][j][k][l]=(-h/(2*m))*Laplacian(temp_z_u,i,j,k,l)+V(i,j,k,l)*temp_z_u[i][j][k][l];

                        J4_v[i][j][k][l]=(h/(2*m))*Laplacian(temp_x_v,i,j,k,l)+V(i,j,k,l)*temp_x_v[i][j][k][l];
                        K4_v[i][j][k][l]=(h/(2*m))*Laplacian(temp_y_v,i,j,k,l)+V(i,j,k,l)*temp_y_v[i][j][k][l];
                        L4_v[i][j][k][l]=(h/(2*m))*Laplacian(temp_z_v,i,j,k,l)+V(i,j,k,l)*temp_z_v[i][j][k][l];
                    }
                }                            
             }
          }



          //Update The Points on the u and v grids repspectively 
          
          for(int i=0; i<2; i++){
            for(int j=0; j<size; j++){
                for(int k=0; k<size; k++){
                    for(int l=0; l<size; l++){
                        u_x[i][j][k][l]+=(h_step/6)*(J1_u[i][j][k][l]+(2*J2_u[i][j][k][l])+(2*J3_u[i][j][k][l])+J4_u[i][j][k][l]);
                        u_y[i][j][k][l]+=(h_step/6)*(K1_u[i][j][k][l]+(2*K2_u[i][j][k][l])+(2*K3_u[i][j][k][l])+K4_u[i][j][k][l]);
                        u_z[i][j][k][l]+=(h_step/6)*(L1_u[i][j][k][l]+(2*L2_u[i][j][k][l])+(2*L3_u[i][j][k][l])+L4_u[i][j][k][l]);

                        v_x[i][j][k][l]+=(h_step/6)*(J1_v[i][j][k][l]+(2*J2_v[i][j][k][l])+(2*J3_v[i][j][k][l])+J4_v[i][j][k][l]);
                        v_y[i][j][k][l]+=(h_step/6)*(K1_v[i][j][k][l]+(2*K2_v[i][j][k][l])+(2*K3_v[i][j][k][l])+K4_v[i][j][k][l]);
                        v_z[i][j][k][l]+=(h_step/6)*(L1_v[i][j][k][l]+(2*L2_v[i][j][k][l])+(2*L3_v[i][j][k][l])+L4_v[i][j][k][l]);
                    }
                }
             }
          }

          for(int i=0; i<2; i++){
            for(int j=0; j<size; j++){
                for(int k=0; k<size; k++){
                    for(int l=0; l<size; l++){
                        System.out.println();
                        System.out.printf("Vector u = <%f,%f,%f> \t",u_x[i][j][k][l],u_y[i][j][k][l],u_z[i][j][k][l]);
                        System.out.printf("Vector v = <%f,%f,%f>",v_x[i][j][k][l],v_y[i][j][k][l],v_z[i][j][k][l]);
                        System.out.println();
                    }
                }
            }
            }

        }
    }

public class ROSEMARY{
    public static void main(String[] args){
        TimeIndependentSchrodinger TS=new TimeIndependentSchrodinger();
        TS.UpdateGrids();
    }
}
