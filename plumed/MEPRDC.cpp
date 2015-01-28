/*
 Collective variable for implementation of Residual Dipolar Coupling restraints according
 to the maximum entropy principle as described in Olsson et al 2015.
 
 requires boost for sampling of random magnetic field orientations. It is possible to use 
 the random number generator implemented in PLUMED2 as well - it requires some modification of
 the code below.

  S OLSSON
  A CAVALLI

*/

#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
//#include "tools/Random.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_on_sphere.hpp>
#include <boost/random/variate_generator.hpp>

using namespace std;

namespace PLMD{

   class r_on_sphere {
      public:
         typedef boost::mt19937 gen_type;
         gen_type rand_gen;
         boost::uniform_on_sphere<double> unif_sphere;
         boost::variate_generator< gen_type*, boost::uniform_on_sphere<double> >    m_random_on_sphere;
         r_on_sphere(int seed):
            unif_sphere(3),m_random_on_sphere(&rand_gen, unif_sphere){
               rand_gen.seed(seed);
            }

         vector<double> point(){return m_random_on_sphere();}
   };

   //+PLUMEDOC COLVAR MEPRDC 
   /*
      This collective variable implements Maximum Entropy RDCs according to Olsson et al 2015
      */
   //+ENDPLUMEDOC

   class MEPRDC : public Colvar {
      unsigned  ens_dim;
      bool ascale;
      float dt;
      float sdmax;
      unsigned nave,naveps;
      bool w_period;
      unsigned count,scount,nadds;
      vector<double> a;
      vector<double> coupl;
      vector<double> rdc_mean;
      vector<double> lambda;
      vector<unsigned> testset;
      string string_prefix;
      public:
      r_on_sphere *ros;
      MEPRDC(const ActionOptions&);
      static void registerKeywords( Keywords& keys );
      virtual void calculate();
      double energy(vector<double>);
      vector<double> uniform_on_sphere();
      void sample_a_conditional(bool upd_);
   };

   PLUMED_REGISTER_ACTION(MEPRDC,"MEPRDC")

   void MEPRDC::registerKeywords( Keywords& keys ){
      Colvar::registerKeywords( keys );
      keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose.");
      keys.add("numbered","ATOMS","the couple of atoms involved in each of the bonds for which you wish to calculate the RDC. "
            "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one dipolar coupling will be "
            "calculated for each ATOMS keyword you specify.");
      keys.reset_style("ATOMS","atoms");
      keys.add("numbered","COUPLING","Add an experimental value for each coupling. ");
      keys.add("numbered","TESTSET","Add a flag to indicate whether or not coupling is in the test-set. (1: in testset, 0: in training set)");
      keys.addFlag("WRITE_RDC",true,"Write the back-calculated RDCs and estimated parameters after every estimation"); 
      keys.add("compulsory","OUTPUT_PREFIX","RDC_OUT","Prefix for output files");
      keys.add("compulsory","SEED","Seed for RNG");
      keys.add("compulsory","DMAXS","50","Initialize D_max*s");
      keys.add("compulsory","DT","0.0001","Set stepsize for parameter estimation");
      keys.add("compulsory","NAVERAGE","1000","Stride for the collection of statistics.");
      keys.add("compulsory","NAVEPS","100","How many strides per estimation.");
      keys.addFlag("AUTOSCALE",true,"whether or not to optimize the degree of alignment parameter, s");
      keys.remove("NOPBC");
   }

   MEPRDC::MEPRDC(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
   {
      vector<AtomNumber> t, atoms;
      for(unsigned i=1;;++i ){
         parseAtomList("ATOMS", i, t );
         if( t.empty() ) break;
         if( t.size()!=2 ){
            string ss; Tools::convert(i,ss);
            error("ATOMS" + ss + " keyword has the wrong number of atoms");
         }
         atoms.push_back(t[0]);
         atoms.push_back(t[1]);
         t.resize(0);
      }
      // Read in RDC values
      coupl.resize( atoms.size()/2 );
      rdc_mean.resize( atoms.size()/2 );
      lambda.resize( atoms.size()/2,0.0 );
      testset.resize( atoms.size()/2, 0 );
      unsigned ntarget=0;
      for(unsigned i=0;i<coupl.size();++i){
         if( !parseNumbered( "COUPLING", i+1, coupl[i] ) ) break;
         ntarget++;
      }
      if( ntarget!=coupl.size() ) error("found wrong number of COUPLING values");
      ntarget=0;
      for(unsigned i=0;i<coupl.size();++i){
         if( !parseNumbered( "TESTSET", i+1, testset[i] ) ) break;
         ntarget++;
      }
      if( ntarget!=coupl.size() ) error("found wrong number of TESTSET values");
      ascale=true;
      dt=0.0001;
      sdmax=50.;
      nave=10;
      naveps=100;
      w_period=false;
      unsigned seed=0;
      parseFlag("WRITE_RDC", w_period);
      parse("SEED", seed);

      string_prefix = "rdc";
      parse("OUTPUT_PREFIX", string_prefix);

      parse("DMAXS",sdmax);
      parse("DT",dt);
      parseFlag("AUTOSCALE", ascale);
      parse("NAVERAGE", nave); 
      parse("NAVEPS", naveps); 

      checkRead();

      log<<"  Bibliography "
         <<plumed.cite("Olsson S; Ekonomiuk D; Sgrignani J; Cavalli A 'Molecular dynamics of biomolecules through direct analysis of dipolar couplings' Submitted. ") <<"\n";

      checkRead();
      addValueWithDerivatives();
      setNotPeriodic();
      requestAtoms(atoms);
      ros = new r_on_sphere( seed+multi_sim_comm.Get_rank() );
      a = ros->point();
      count = 0;
      scount = 0;
      nadds = 0;
      for(unsigned i=0;i<coupl.size();++i){
         log.printf(" The %dth Bond Dipolar Coupling is calculated from atoms : %d %d.", i+1, atoms[2*i].serial(), atoms[2*i+1].serial());
         log.printf(" Dipolar Coupling is %f.\n",coupl[i]);
      }
      log.printf("  DONE!\n"); log.flush();
   }


   void MEPRDC::sample_a_conditional(bool upd_){
      double E = energy(a);
      vector<double> a_buff = a;
      for(unsigned i=0;i<10;i++){
         //a = uniform_on_sphere();
         a = ros->point();
         double Etmp  = energy(a);

         bool accept = false;

         if(Etmp<=E){
            accept = true;
         }
         else {
            double r = ros->rand_gen();
            if(r<std::exp((E-Etmp))){
               accept = true;
            }
         }

         if(accept){
            E = Etmp;
            a_buff = a;
         }
         else
            a = a_buff;
      }

      if(upd_){

         int mpi_rank = multi_sim_comm.Get_rank();
         int mpi_size = multi_sim_comm.Get_size();


         //SEND
         double rdc[rdc_mean.size()];
         double lambda_mpi[lambda.size()+1];
         double nrep = 1./(double)mpi_size;
         for(unsigned i=0;i<rdc_mean.size();i++){
            rdc[i] = rdc_mean[i]/(double)nadds;
         }
         if(comm.Get_rank()==0) {
            multi_sim_comm.Sum(&rdc[0], coupl.size() );
            for(unsigned i=0;i<coupl.size();i++) rdc[i] *= nrep;
         } else for(unsigned i=0;i<coupl.size();i++) rdc[i] = 0.;
         comm.Sum(&rdc[0], coupl.size());

         if(mpi_rank==0) {
            for(int i=0;i<rdc_mean.size();i++) {
               rdc_mean[i] = rdc[i];
            }

            double s = 1;
            if(ascale){
               //scaling mean
               double D2 = 0;
               double DD = 0;
               for(unsigned i=0;i<lambda.size();i++){
                  if(testset[i]!=1) {
                     D2 += sdmax*sdmax*rdc_mean[i]*rdc_mean[i];
                     DD += sdmax*rdc_mean[i]*coupl[i];
                  }
               }
               s = DD/D2;
               if(s>1.2) s = 1.2;
               if(s<0.8) s = 0.8;
               sdmax = s*sdmax;
            }


            bool printout=false;
            FILE *outfile, *outfile2;
            if(w_period) printout = (!(getStep()%w_period));
            if(printout) {
               outfile = fopen((string_prefix+string("_rdc.val")).c_str(), "w");
               for(unsigned r=0;r<coupl.size();r++) {
                  fprintf(outfile," %4i %10.6lf %10.6lf\n", r+1, sdmax*rdc_mean[r], coupl[r]);
               }
               fclose(outfile);

               outfile2 = fopen((string_prefix+string("_rdc.dat")).c_str(), "a");
               for(unsigned r=0;r<coupl.size();r++) {
                  fprintf(outfile2,"%10.6lf ", sdmax*rdc_mean[r]);
               }
               fprintf(outfile2,"\n");
               fclose(outfile2);
            }

            if(printout) {
               outfile = fopen((string_prefix+string("_Lambda.val")).c_str(), "w");
               outfile2 = fopen((string_prefix+string("_Lambda.dat")).c_str(), "a");
            }

            double ss = 0, rr = 0;
            double sstest = 0, rrtest = 0;
            bool test_set_used=false;
            for(unsigned i=0;i<lambda.size();i++){
               if(testset[i]!=1) {
                  // ! estimate lambda if not in testset
                  lambda[i] = lambda[i]+dt*(sdmax*rdc_mean[i] - coupl[i]);
                  ss += (sdmax*rdc_mean[i] - coupl[i])*(sdmax*rdc_mean[i] - coupl[i]);
                  rr += coupl[i]*coupl[i];
               } else {
                  sstest += (sdmax*rdc_mean[i] - coupl[i])*(sdmax*rdc_mean[i] - coupl[i]);
                  rrtest += coupl[i]*coupl[i];
                  test_set_used=true;
               }
               if(!test_set_used) {
                  sstest=-1.; rrtest=1.;
               }
               lambda_mpi[i] = lambda[i];
               if(printout) { 
                  fprintf(outfile2,"%10.6lf ", sdmax*lambda[i]);
                  fprintf(outfile,"%10.6lf\n ", lambda[i]);
               }
            }
            lambda_mpi[lambda.size()] = sdmax;
            if(printout) {
               fprintf(outfile2,"\n");
               fclose(outfile);
               fclose(outfile2);
            }
            if(ascale){
               log.printf("Updating lambdas: [ %f ] [ %f ][ %f ][ %f ]\n", sqrt(ss/rr), sqrt(sstest/rrtest), s, sdmax); 
            }
            else {
               log.printf("Updating lambdas: [ %f ] [ %f ]\n", sqrt(ss/rr), sqrt(sstest/rrtest)); 
            }
         }
         multi_sim_comm.Bcast(&lambda_mpi[0], lambda.size()+1, 0);
         //comm.Sum(&lambda_mpi[0], coupl.size()+1);
         for(unsigned i=0;i<lambda.size();i++){
            lambda[i] = lambda_mpi[i];
            rdc_mean[i] = 0;
         }
         nadds = 0;
         if(ascale && mpi_rank!=0) 
            sdmax = lambda_mpi[lambda.size()];
         log.flush();
      }
   }


   double MEPRDC::energy(vector<double> abuff) {
      unsigned N = getNumberOfAtoms();
      double E = 0;
      unsigned index=0; 
      for(unsigned i=0;i<N;i+=2){ // loop over atom-pairs 
         Vector dxdydz;
         dxdydz = pbcDistance(getPosition(i),getPosition(i+1));
         double d = 1./dxdydz.modulo();
         double cosa = dxdydz[0]*abuff[0] + dxdydz[1]*abuff[1] + dxdydz[2]*abuff[2];
         cosa = cosa*d;
         double p2 = 0.5*(3.0*cosa*cosa-1.0);
         E = E + sdmax*p2*lambda[index];
         ++index;
      }
      return E;
   }

   void MEPRDC::calculate()
   {
      unsigned N = getNumberOfAtoms();
      vector<Vector> deriv(N);
      double E = 0;
      E = 0;
      bool add_=false;
      bool upd_=false;
      ++count;
      ++scount;

      if(count==nave*naveps) {
         count=0;
         upd_=true;
      }
      if(scount==nave) {
         scount = 0;
         add_=true;
         ++nadds;
      }
      unsigned index=0; 
      for(unsigned i=0;i<N;i+=2){ // loop over atom-pairs 
         Vector dxdydz;
         dxdydz = pbcDistance(getPosition(i),getPosition(i+1));
         double d = dxdydz.modulo();
         double invd = 1./d;
         double dx = dxdydz[0];
         double dy = dxdydz[1];
         double dz = dxdydz[2];
         double cosa  = dx*a[0] + dy*a[1] + dz*a[2];
         cosa = cosa*invd;
         double p2 = 0.5*(3.0*cosa*cosa-1.0);
         E = E + sdmax*p2*lambda[index];
         double fact = cosa*invd;         
         double tmp = sdmax*lambda[index]*3.0*cosa;
         deriv[i][0] = tmp*(a[0]*invd - fact*dx);
         deriv[i][1] = tmp*(a[1]*invd - fact*dy);
         deriv[i][2] = tmp*(a[2]*invd - fact*dz);

         deriv[i+1] = -deriv[i];
         if(add_) {
            rdc_mean[index] =  rdc_mean[index] + p2;
         }
         ++index;
      }
      sample_a_conditional(upd_);
      Tensor virial;
      virial.zero();
      Vector For;
      for(unsigned i=0;i<N;i++)
      {
         For=deriv[i];
         setAtomsDerivatives(i,For);
         virial=virial+(-1.*Tensor(getPosition(i),For));
      }

      setValue           (E);
      setBoxDerivatives  (virial);

   }

}
