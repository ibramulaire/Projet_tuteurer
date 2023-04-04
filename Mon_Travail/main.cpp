/*!
 * \file main.cpp
 * \brief Ce Programme permet de decimer un maillage en utilisant la saillance .
Pour compiler, il faut dans un premier temps instaler les librairis pmp https://www.pmp-library.org/ ,assimp et armadillo.

Pour l'utilisation, il faut passer en parametre le fichier à traité, le nombre de contraction d'arrete à faire et un nombre entre 0 et 100.
Ce dernier permet de specifier le poid de la saillance.
Par exemple:
Si ce nombre est de 50 alors la saillance aura le meme poid que l'erreur quadratique. 
S'il est de 20, l'ereur quadratique aura alors 80 de poid.  

 * \author \Diarra Ibrahim
 * \version 0.1
 */
#include <iostream>
#include <string>
#include <vector>
#include <math.h> 
#include<pmp/io/io.h>
#include <pmp/SurfaceMesh.h>
#include <pmp/Types.h>
#include <pmp/algorithms/Curvature.h>
#include<pmp/BoundingBox.h>
#include<pmp/utilities.h>
#include <pmp/algorithms/Normals.h>
#include <pmp/MatVec.h>
#include "Quadric.h"
#include<pmp/algorithms/Quadric.h>
using namespace std;

/*!
 * \fn double  compute_diagonalLength(const pmp::SurfaceMesh & ppmmesh)
 * \brief Cette fonction calcule la longeur de la diagonale du maillage trait.
 *
 * \param ppmmesh Reference un object SurfaceMesh de la biliotheque pmp(polynome mesh processing) qui est un maillage .
 * \return retourne la valeur calculer dans un double.
 */
double  compute_diagonalLength(const pmp::SurfaceMesh & ppmmesh)
{
     double xmin ,xmax,ymin,ymax,zmin ,zmax;
     xmin=ymin=zmin=std::numeric_limits<double>::max();
     xmax=ymax=zmax=std::numeric_limits<double>::min();
     pmp::VertexProperty<pmp::Point>  points = ppmmesh.get_vertex_property<pmp::Point>("v:point");

     for (auto v : ppmmesh.vertices())     
     {
          auto p= points[v];
          if(p[0]<xmin)
          xmin=p[0];
          if(p[1]<ymin)
          ymin=p[1];
          if(p[2]<zmin)
          zmin=p[2];

          if(p[0]>xmax)
          xmax=p[0];
          if(p[1]>ymax)
          ymax=p[1];
          if(p[2]>zmax)
          zmax=p[2];
          
     }

     return sqrt(pow(xmax-xmin,2)+pow(ymax-ymin,2)+pow(zmax-zmin,2));

}
/*!
 * \fn void computeSaillancy(pmp::SurfaceMesh & ppmmesh)
 * 
 * \brief Cette fonction calcule la saillance en chaque sommet du maillage traité avce la .
 *
 * \param ppmmesh Reference un object SurfaceMesh de la biliotheque pmp(polynome mesh processing) qui est un maillage .
 * \return void
 */
void computeSaillancy(pmp::SurfaceMesh & ppmmesh)
{
     pmp::Curvature c= pmp::Curvature(ppmmesh);
     c.analyze_tensor();

     
     pmp::VertexProperty<pmp::Point>  points = ppmmesh.get_vertex_property<pmp::Point>("v:point");

     auto saillancyi = ppmmesh.add_vertex_property<vector<pmp::Scalar>>("added:saillancyi");

     auto smoothsaillancy = ppmmesh.add_vertex_property<pmp::Scalar>("added:smoothsaillancy");

     vector<double> max_saillancyi;

     for(int i=0;i<5;i++)
          max_saillancyi.push_back(std::numeric_limits<double>::min());

    

     double diagonalLength=compute_diagonalLength(ppmmesh);


     double epsilon=2*0.003*diagonalLength;
     vector<int> indicesepsilon={2,3,4,5,6,8,10,12};
     for (auto v : ppmmesh.vertices())     
     {
       
          vector<pmp::Vertex> voisinnages2;
          vector<pmp::Vertex> voisinnages3;
          vector<pmp::Vertex> voisinnages4;
          vector<pmp::Vertex> voisinnages5;
          vector<pmp::Vertex> voisinnages6;
          vector<pmp::Vertex> voisinnages8;
          vector<pmp::Vertex> voisinnages10;
          vector<pmp::Vertex> voisinnages12;
          for(auto vn: ppmmesh.vertices())
          {
               
               double distance=sqrt((pow(points[vn][0]-points[v][0],2)+pow(points[vn][1]-points[v][1],2)+pow(points[vn][2]-points[v][2],2)));
               if(distance<2*epsilon)
                    voisinnages2.push_back(vn);
               if(distance<3*epsilon)
                    voisinnages3.push_back(vn);  
               if(distance<4*epsilon)
                    voisinnages4.push_back(vn);
               if(distance<5*epsilon)
                    voisinnages5.push_back(vn); 
               if(distance<6*epsilon)
                    voisinnages6.push_back(vn);
               if(distance<8*epsilon)
                    voisinnages8.push_back(vn);  
               if(distance<10*epsilon)
                    voisinnages10.push_back(vn);
               if(distance<12*epsilon)
                    voisinnages12.push_back(vn); 
          }
          vector<vector<pmp::Vertex>> finalv;
          finalv.push_back(voisinnages2);
          finalv.push_back(voisinnages3);
          finalv.push_back(voisinnages4);
          finalv.push_back(voisinnages5);
          finalv.push_back(voisinnages6);
          finalv.push_back(voisinnages8);
          finalv.push_back(voisinnages10);
          finalv.push_back(voisinnages12);
       
 
          epsilon=0.003*diagonalLength;
          vector<pmp::Scalar> gaussianWeigths;
          for(int i=0;i<8;i++)  
          {
          
          double numerateur=0;   
          double denominateur =0;

          double sigma=pow(indicesepsilon[i]*epsilon,2);
          
          for(int j=0;j<(int)finalv[i].size();j++)
          {
               double norme=0;
               double expo=0;
               norme=(pow(points[finalv[i][j]][0]-points[v][0],2)+pow(points[finalv[i][j]][1]-points[v][1],2)+pow(points[finalv[i][j]][2]-points[v][2],2));
               expo=exp(-norme/2*sigma);
               numerateur+=c.mean_curvature(finalv[i][j])*expo;
               denominateur+=expo;
               
          }
          
          if(denominateur==0||numerateur==0)
               gaussianWeigths.push_back(0.0);
          else
               gaussianWeigths.push_back((numerateur/denominateur));

          }
           vector<pmp::Scalar> saillancy;
           double s1,s2,s3,s4,s5;
                                 
          s1=abs(gaussianWeigths[0]-gaussianWeigths[2]);
          s2=abs(gaussianWeigths[1]-gaussianWeigths[4]);
          s3=abs(gaussianWeigths[2]-gaussianWeigths[5]);
          s4=abs(gaussianWeigths[3]-gaussianWeigths[6]);
          s5=abs(gaussianWeigths[4]-gaussianWeigths[7]);
          
          max_saillancyi[0]=max( max_saillancyi[0],s1);
          max_saillancyi[1]=max( max_saillancyi[1],s2);
          max_saillancyi[2]=max( max_saillancyi[2],s3);
          max_saillancyi[3]=max( max_saillancyi[3],s4);
          max_saillancyi[4]=max( max_saillancyi[4],s5);

           saillancy.push_back(s1);
           saillancy.push_back(s2);
           saillancy.push_back(s3);
           saillancy.push_back(s4);
           saillancy.push_back(s5);
           
       
          saillancyi[v]= saillancy;   
           
          }
          
          double maxSmmothsaillancy=0; 

     
     for (auto v : ppmmesh.vertices())     
     {
	     vector<double> local_max_saillancyi;
	     for(int i=0;i<5;i++)
		  local_max_saillancyi.push_back( saillancyi[v][i]);
	     
	     for (auto vn : ppmmesh.vertices(v))     
	     {
		  for(int i=0;i<5;i++){
		  local_max_saillancyi[i] = max(local_max_saillancyi[i], (double)saillancyi[vn][i]);
		  
		  }
	     }
	     
	     double finalsaillancy=0;
	     double saliencySum=0;

	     for(int i=0;i<5;i++)
	     {
		  double   factor = pow(max_saillancyi[i]-local_max_saillancyi[i],2);
		 
		  finalsaillancy += saillancyi[v][i] * factor;
		  saliencySum += factor;
	     }
	     if(finalsaillancy==0||saliencySum==0)
	      smoothsaillancy[v] =0;
	     else
	     {
	      smoothsaillancy[v] =finalsaillancy/saliencySum;
	      if(maxSmmothsaillancy<smoothsaillancy[v] )
	         maxSmmothsaillancy= smoothsaillancy[v] ;
	    // cout<< smoothsaillancy[v]<<endl;
	     }
     }
     
     for (auto v : ppmmesh.vertices())     
     {
      smoothsaillancy[v]= smoothsaillancy[v]/maxSmmothsaillancy;
     }
     
}

/*!
 * \fn void computeQuadrix(pmp::SurfaceMesh & ppmmesh)
 * 
 * \brief Cette fonction sert à la simplification. Elle permet de calculer les quadrics de depart.
 *
 * \param ppmmesh Reference un object SurfaceMesh de la biliotheque pmp(polynome mesh processing) qui est un maillage .
 * \return void
 */
void computeQuadrix(pmp::SurfaceMesh & ppmmesh)
{
  pmp::Normals::compute_face_normals(ppmmesh);
 // auto quadrics =ppmmesh.add_vertex_property<pmp::Scalar>("added:quadrics");
  auto normales= ppmmesh.get_face_property<pmp::Normal>("f:normal");
  pmp::VertexProperty<pmp::Point>  points = ppmmesh.get_vertex_property<pmp::Point>("v:point");
  auto Qmat =ppmmesh.add_vertex_property<arma::mat>("added:Qmat");

 for (auto v:ppmmesh.vertices()) 
  {
      
      Quadric q({0,0,0},{0,0,0});
      for(auto f:ppmmesh.faces(v))
      {
        // cout<<Quadric({normales[f][0],normales[f][1],normales[f][2]},{points[v][0],points[v][1],points[v][1]}).error({points[v][0],points[v][1],points[v][1],1})<<endl;
       q.add(Quadric({normales[f][0],normales[f][1],normales[f][2]},{points[v][0],points[v][1],points[v][2]}));
      }
      Qmat[v]=q.getMat();
      
    // cout<<q.error({points[v][0],points[v][1],points[v][2],1})<<endl;
  }
  
}

/*!
 * \fn arma::vec Optimalposition(arma::mat qprime)
 * 
 * \brief Cette fonction sert à la simplification. Elle calcule la position optimal en minimisant la fonction d'erreur .
 *
 * \param ppmmesh Reference un object SurfaceMesh de la biliotheque pmp(polynome mesh processing) qui est un maillage .
 * \return return un le point optimal
 */
arma::vec Optimalposition(arma::mat qprime)
{
      arma::mat res= qprime;
      //cout<<res<<endl;
      res(3,0)=0;
      res(3,1)=0;
      res(3,2)=0;
      res(3,3)=1;
     // cout<<res<<endl;
      res= arma::inv(res);
  //    cout<<res<<endl;
     return{res(0,3),res(1,3),res(2,3),res(3,3)};


}
/*!
 * \fn void ComputePAirInfos(pmp::SurfaceMesh & ppmmesh,int poid)
 * 
 * \brief Cette fonction sert à la simplification. Elle calcule la position optimal en minimisant la fonction d'erreur .
 *
 * \param ppmmesh Reference un object SurfaceMesh de la biliotheque pmp(polynome mesh processing) qui est un maillage .
 * 
 * \param poid Poid de la saillance par rapport à l'erreur quadratique moyenne .
 */
void ComputePAirInfos(pmp::SurfaceMesh & ppmmesh,int poid)
{
     cout<< "Compute pair infos"<<endl;
     auto quadrics =ppmmesh.get_vertex_property<arma::mat>("added:Qmat");
     auto Qprime =ppmmesh.add_edge_property<arma::mat>("added:Qprime");
     auto Optimaleposition =ppmmesh.add_edge_property<arma::vec>("added:optimal");
     auto Cout =ppmmesh.add_edge_property<pmp::Scalar>("added:cout");
     auto Cout2 =ppmmesh.add_edge_property<pmp::Scalar>("added:cout2");
     auto Cout3 =ppmmesh.add_edge_property<pmp::Scalar>("added:cout3");
     auto smoothsaillancy = ppmmesh.get_vertex_property<pmp::Scalar>("added:smoothsaillancy");
     for (auto e:ppmmesh.edges()) 
     {
          Qprime[e]=quadrics[ppmmesh.vertex(e,0)]+quadrics[ppmmesh.vertex(e,1)];
          Optimaleposition[e]=Optimalposition( Qprime[e]);
          Quadric q(Qprime[e]);
          Cout[e]=q.error(Optimaleposition[e]);
         // cout<< Cout[e]<<endl;
     }


}
/*!
 * \fn  pmp::Edge bestCout(pmp::SurfaceMesh & ppmmesh,int poid)
 * 
 * \brief Cette fonction sert à la simplification. 
 * Elle calcule l'arete avec le coût de decimation le plus faible en fonction de la saillance et de l'erreur quadratique  .
 *
 * \param ppmmesh Reference un object SurfaceMesh de la biliotheque pmp(polynome mesh processing) qui est un maillage .
 * 
 * \param poid Poid de la saillance par rapport à l'erreur quadratique moyenne .
 */
 pmp::Edge bestCout(pmp::SurfaceMesh & ppmmesh,int poid)
 {
     
      auto Cout =ppmmesh.get_edge_property<pmp::Scalar>("added:cout");
      auto Cout2 =ppmmesh.get_edge_property<pmp::Scalar>("added:cout2");
      auto Cout3 =ppmmesh.get_edge_property<pmp::Scalar>("added:cout3");
      auto smoothsaillancy = ppmmesh.get_vertex_property<pmp::Scalar>("added:smoothsaillancy");
      double maxcout=0;
      for (auto e:ppmmesh.edges()) 
          {
               if(Cout[e]>maxcout)
               maxcout=Cout[e];
          }

      for (auto e:ppmmesh.edges()) 
          {
               Cout3[e] =Cout[e]/maxcout;
               Cout2[e]=((100-poid)*(Cout[e]/maxcout)+1)*(poid*((smoothsaillancy[ppmmesh.vertex(e,0)]+smoothsaillancy[ppmmesh.vertex(e,1)]))+1);
          }
     

      
      double bestcost=std::numeric_limits<double>::max();
      pmp::Edge best;
     for (auto e:ppmmesh.edges()) 
          {
           if( ppmmesh.is_collapse_ok(ppmmesh.find_halfedge(ppmmesh.vertex(e,0),ppmmesh.vertex(e,1) )))
           {
               
               if(Cout2[e]==0)
               return e;
               else
               if(Cout2[e]<bestcost)
               {
                    bestcost=Cout2[e];
                    best=e;
               }
           }
          
          }
     return best;

 }


/*!
 * \fn void updatePaireInfo(pmp::SurfaceMesh & ppmmesh,pmp::Vertex & toupdate,int poid)
 * 
 * \brief Cette fonction sert à la simplification. Elle permet de calculer de mettre a jour les quadrics apres une contraction.
 *
 * \param ppmmesh Reference un object SurfaceMesh de la biliotheque pmp(polynome mesh processing) qui est un maillage .
 * \param poid Poid de la saillance par rapport à l'erreur quadratique moyenne .
 * \return void
 */
 void updatePaireInfo(pmp::SurfaceMesh & ppmmesh,pmp::Vertex & toupdate,int poid)
 { 
     pmp::Normals::compute_face_normals(ppmmesh);
   //auto quadrics =ppmmesh.add_vertex_property<pmp::Scalar>("added:quadrics");
     auto normales= ppmmesh.get_face_property<pmp::Normal>("f:normal");
     pmp::VertexProperty<pmp::Point>  points = ppmmesh.get_vertex_property<pmp::Point>("v:point");
     auto Qmat =ppmmesh.get_vertex_property<arma::mat>("added:Qmat");
     
     auto Qprime =ppmmesh.get_edge_property<arma::mat>("added:Qprime");
     auto Optimaleposition =ppmmesh.get_edge_property<arma::vec>("added:optimal");
     auto Cout =ppmmesh.get_edge_property<pmp::Scalar>("added:cout");
     auto smoothsaillancy = ppmmesh.get_vertex_property<pmp::Scalar>("added:smoothsaillancy");
   
     for (auto v:ppmmesh.vertices(toupdate)) 
     {

      Quadric q({0,0,0},{0,0,0});
      for(auto f:ppmmesh.faces(v))
      {
   
       q.add(Quadric({normales[f][0],normales[f][1],normales[f][2]},{points[v][0],points[v][1],points[v][2]}));
      }
      Qmat[v]=q.getMat();

      pmp::Edge e=ppmmesh.find_edge(toupdate,v);

      Qprime[e]=Qmat[ppmmesh.vertex(e,0)]+Qmat[ppmmesh.vertex(e,1)];
  
      Optimaleposition[e]=Optimalposition( Qprime[e]);
  
      Quadric Q(Qprime[e]);
 
      Cout[e]=Q.error(Optimaleposition[e]);
     }

 }




/*!
 * \fn void decimer(pmp::SurfaceMesh & ppmmesh,int poid)
 * 
 * \brief Cette fonction sert à la simplification. Elle permet de contracter l'arrete avec le coût le moins elever .
 *
 * \param ppmmesh Reference un object SurfaceMesh de la biliotheque pmp(polynome mesh processing) qui est un maillage .
 * \param poid Poid de la saillance par rapport à l'erreur quadratique moyenne .
 * \return void
 */

void decimer(pmp::SurfaceMesh & ppmmesh,int poid)
{
     
     auto quadrics =ppmmesh.get_vertex_property<arma::mat>("added:Qmat");
     auto Qprime =ppmmesh.get_edge_property<arma::mat>("added:Qprime");
     auto Optimaleposition =ppmmesh.get_edge_property<arma::vec>("added:optimal");
     auto Cout =ppmmesh.get_edge_property<pmp::Scalar>("added:cout");
     auto Cout2 =ppmmesh.get_edge_property<pmp::Scalar>("added:cout2");
     auto Cout3 =ppmmesh.get_edge_property<pmp::Scalar>("added:cout3");
     auto smoothsaillancy = ppmmesh.get_vertex_property<pmp::Scalar>("added:smoothsaillancy");
     pmp::Edge e= bestCout(ppmmesh,poid);
     pmp::Vertex v0=ppmmesh.vertex(e,0);
     pmp::Vertex v1=ppmmesh.vertex(e,1);
     quadrics[v1]=  Qprime[e];
     smoothsaillancy[v1]=(smoothsaillancy[v0]+smoothsaillancy[v1]);
   //  cout<<"cout:"<<Cout[e]<<"cout2:"<<Cout2[e]<<"cout3:"<<Cout3[e]<<" saillancy:"<< smoothsaillancy[v1]<<endl;
     ppmmesh.position(v1)[0]= Optimaleposition[e](0);
     ppmmesh.position(v1)[1]= Optimaleposition[e](1);
     ppmmesh.position(v1)[2]= Optimaleposition[e](2);
     ppmmesh.position(v0)[0]= Optimaleposition[e](0);
     ppmmesh.position(v0)[1]= Optimaleposition[e](1);
     ppmmesh.position(v0)[2]= Optimaleposition[e](2);
     ppmmesh.collapse( ppmmesh.find_halfedge(v0,v1 ));
   

updatePaireInfo(ppmmesh,v1,poid);

}



int main(int argc, char const *argv[])
{
      string filename;
      string destdir;
      string meshdir;
      int nb;
      int poid;

     filename =string(argv[1]); 
     meshdir = string(argv[2]); 
     destdir = string(argv[5]);
     nb=stoi(string(argv[3]));
     poid=stoi(string(argv[4]));
    
     pmp::SurfaceMesh ppmmesh;
     pmp::read(ppmmesh,meshdir+filename);
    
    
     if(nb>99 || nb<1)
     {
          cout<<"le troisieme parametre doit un nombre entre 2 et 99 et represente le pourcentatge de decimation souhaité "<<endl;
            return 0;
     }
     
      nb=(ppmmesh.n_vertices()* nb/100);
      if( nb>=ppmmesh.n_vertices())
      nb--;
     
     if(poid>100)
     {
          cout<<"le quatrieme parametre doit etre compris entre 1 et 10"<<endl;
            return 0;
     }
     computeSaillancy(ppmmesh);

     computeQuadrix(ppmmesh);
     ComputePAirInfos(ppmmesh,poid);

     for(int i=1;i<=nb;i++)
     {
      decimer( ppmmesh,poid);
      cout<<i<<" "<<(double)i*100/(double)nb<<"%"<<endl;
     }
    ppmmesh.garbage_collection();




    
 
   string outputname=destdir+string("/nbD:")+string(argv[3])+string("SPoi")+string(argv[4])+filename;
    
   //  auto smoothsaillancy = ppmmesh.get_vertex_property<pmp::Scalar>("added:smoothsaillancy");
    // auto quadrics =ppmmesh.get_vertex_property<pmp::Scalar>("added:quadrics");
     

     pmp::write(ppmmesh,outputname);


     return 0;
}


