#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Q2baranov(NumericVector e, List inputs, IntegerVector idx) {
  
  // fleet inputs
  List fleet_inputs = inputs["fleet_inputs"];
  NumericVector fixedEffort = fleet_inputs["effort"];
  int n_fl_mt = fleet_inputs["n_fl_mt"];
  NumericVector eshare = fleet_inputs["eshare"];
  IntegerVector fl_idx = fleet_inputs["fl_idx"];
  
  
  // combine estimated and fixed effort in a single vector and multiply with effort shares 
  NumericVector E(n_fl_mt); 
  int ctr_e = 0; 
  int fl_idx_i_1 = 0;
  for(int i = 0; i < n_fl_mt; i++){
    if(idx(fl_idx(i)) >= 0){
      if(fl_idx(i) > fl_idx_i_1){
        ctr_e += 1;
      }
      E(i) = e(ctr_e) * eshare(i);
      fl_idx_i_1 = fl_idx(i);
    } else {
      E(i) = fixedEffort(fl_idx(i)) *  eshare(i);
    }
  } 
  
  List stock_inputs = inputs["stock_inputs"];
  int n_fleets = idx.length();
  int n_e = e.length();
  NumericVector y(n_e);
  
  // compute quota minus catch for each fleet for which effort has to be estimated
  int ctr = 0;
  for(int i = 0; i < n_fleets; i++){
    int id = idx(i);
    // if fleet is not fixedEffort
    if(id >= 0){
      
      // extract info
      List stock = stock_inputs[id];
      NumericVector Q = stock["Q"];    // quota by fleet
      NumericMatrix wt = stock["wt"];  // catch weights-at-age (fleet, age)
      NumericMatrix q = stock["q"];    // catchability-at-age (fleet, age)
      NumericMatrix q_sel = stock["q.sel"];    // catchability-at-age of retained catch part (fleet, age)
      NumericVector M = stock["M"];    // natural mortality rate at age
      NumericVector N = stock["N"];    // stock numbers at age
      
      // fill in the negative value of the quota for fleet i for stock idx[i]
      y(ctr) = -Q(i);
      
      // find métiers
      LogicalVector metier_idx = fl_idx == i ;
      
      // loop over the ages and add the catch for the stock
      int n_ages = N.length();
      double F = 0;
      NumericVector q_e(n_fl_mt);
      for(int j = 0;  j < n_ages; j++){
        // compute F
        F = 0;
        NumericVector q_a = q(_,j);   // vector with catchabilities for all fleet and métiers for age a
        q_e = q_a * E;                // compute partial F by fleet x métier for age a
        F = sum(q_e);     
        
        // métier loop
        for(int k = 0;  k < n_fl_mt; k++){
          if(metier_idx(k)){
            y(ctr) += (wt(k,j) * E(k)*q_sel(k,j)/(F + M(j)) * N(j) * (1-exp(-F - M(j))));  // retained catchabilities
          }
        }
      }
      ctr = ctr + 1;
    }
  }
  return y;
}

// little helper to transform array indices (dim = 6) to vector (one dimensional) index
int flat_index(int dim1, int dim2,int dim3, int dim4, int dim5, int dim6,
               int idx1, int idx2, int idx3, int idx4, int idx5, int idx6){
  return idx1 + idx2 * dim1 + idx3 * dim1 * dim2 +  idx4 * dim1 * dim2 * dim3 +  idx5 * dim1 * dim2 * dim3 * dim4  +  idx6 * dim1 * dim2 * dim3 * dim4 * dim5;
}


// [[Rcpp::export]]
S4 update_catch_effort_FLFleets(NumericMatrix effort,  // new effort matrix [nit x fleets] 
                                S4 fleets,             // FLFleetExts object to be filled
                                List inputs,           // List with for each iteration the stock numbers, landing and discard weights, selectivity, and catcability by stock, fleet and métier
                                int year){             // which year to update
  
  
  S4 cl_fleets = clone(fleets); // make a deep copy
  
  // nr of iterations and fleets
  int nit = effort.nrow();
  int n_fleets = effort.ncol();
  double F = 0;
  
  List inputs0 = inputs[0];
  List fleet_inputs0 = inputs0["fleet_inputs"];
  int n_fl_mt = fleet_inputs0["n_fl_mt"];
  IntegerVector fl_idx = fleet_inputs0["fl_idx"];
  
  
  // extract the data from the fleets object
  List flts = cl_fleets.slot(".Data");
  
  int m_ctr = 0;       // métier counter

  for(int f = 0; f < n_fleets; f++){
    
    // extract single fleet and métiers
    S4 fleet = flts[f];
    S4 metiers = fleet.slot("metiers");
    List mts = metiers.slot(".Data");
    int n_metiers = mts.size();

    S4 flq_effort = fleet.slot("effort");
    
    // convert to vector, fill elements
    IntegerVector dim = flq_effort.attr("dim"); 
    NumericVector effort_dat = flq_effort.slot(".Data");
    
    for(int i = 0; i < nit; i++){
      
      int idx_i = flat_index(1,dim(1),1,1,1,nit, 
                             0, year, 0, 0, 0, i);
      effort_dat(idx_i) = effort(i,f);
      effort_dat.attr("dim") = dim;
      flq_effort.slot(".Data") = effort_dat;
      fleet.slot("effort") = flq_effort;
    }
    
    // loop over métiers
    for(int m = 0; m < n_metiers; m++){
      // extract catch information of each métier
      S4 metier = mts[m];
      S4 catches = metier.slot("catches");
      CharacterVector cnames = catches.slot("names");
      int n_catches = cnames.size();
      List ctchs = catches.slot(".Data");
      
      // loop over catches
      for(int c = 0; c < n_catches;  c++){
        String cname = cnames(c);
        
        // extract landings and discards of single catch object
        S4 catc = ctchs[c];
        
        // Extract objects
        // landings n
        S4 flq_landings_n = catc.slot("landings.n");
        NumericVector landings_n = flq_landings_n.slot(".Data");
        // discard n
        S4 flq_discards_n = catc.slot("discards.n");
        NumericVector discards_n = flq_discards_n.slot(".Data");
        // landings
        S4 flq_landings = catc.slot("landings");
        NumericVector landings = flq_landings.slot(".Data");
        // discards
        S4 flq_discards = catc.slot("discards");
        NumericVector discards = flq_discards.slot(".Data");
        
        // fill across all iterations
        for(int i = 0; i < nit; i++){
          
          List input_i = inputs[i];
          List input_fleet = input_i["fleet_inputs"]; 
          NumericVector eshare = input_fleet["eshare"];
          
          // Vector of efforts by métier
          NumericVector E(n_fl_mt); 
          for(int j = 0; j < n_fl_mt; j++){
            E(j) = effort(i,fl_idx(j)) * eshare(j);
          }
          
          List input_stock = input_i["stock_inputs"]; // list with stock inputs: "N"     "M"     "Q"     "q"     "cw"    "lw"    "dw"    "l.sel"
          
          List stk = input_stock[cname];

          // fill objects
          NumericVector N = stk["N"];
          NumericVector M = stk["M"];
          NumericMatrix q = stk["q"];
          NumericMatrix sel = stk["l.sel"];
          NumericMatrix lw = stk["lw"];
          NumericMatrix dw = stk["dw"];
          
          int n_ages = N.size();
          int idx_i = flat_index(1,dim(1),1,1,1,nit, 
                                 0, year, 0, 0, 0, i);
          double cn = 0;
          double lan = 0;
          double dis = 0;
          NumericVector q_e(n_fl_mt);
          for(int a = 0; a < n_ages; a++){
            F = 0;
            NumericVector q_a = q(_,a);   // vector with catchabilities for all fleet and métiers for age a
            q_e = q_a * E;                // compute partial F by fleet x métier for age a
            F = sum(q_e);                 // compute total F for age a
            
            int idx_ai = flat_index(n_ages,dim(1),1,1,1,nit, 
                                    a, year, 0, 0, 0, i);
            
            cn = E(m_ctr) *q(m_ctr,a)/(F + M(a)) * N(a) * (1-exp(-F - M(a)));
            
            //double cn = 0;
            landings_n[idx_ai] = cn * sel(m_ctr,a);
            discards_n[idx_ai] = cn * (1- sel(m_ctr,a));
            lan +=  (cn * sel(m_ctr,a) * lw(m_ctr,a));
            dis+= (cn * (1- sel(m_ctr,a)) * dw(m_ctr,a)); 
          }
          landings[idx_i] = lan;
          discards[idx_i] = dis; 
        }
        
        // update flquants
        flq_landings_n.slot(".Data") = landings_n;
        flq_landings.slot(".Data") = landings;
        flq_discards_n.slot(".Data") = discards_n;
        flq_discards.slot(".Data") = discards;
        
        catc.slot("landings.n") = flq_landings_n;
        catc.slot("landings") = flq_landings;
        catc.slot("discards.n") = flq_discards_n;
        catc.slot("discards") = flq_discards;
        
        ctchs[c] = catc;
      } // end catch loop
      
      // update objects
      catches.slot(".Data") = ctchs;
      metier.slot("catches") = catches;
      mts[m] = metier;
      
      // next métier
      m_ctr = m_ctr + 1;
    }
    metiers.slot(".Data") = mts;
    fleet.slot("metiers") = metiers;
    
    flts[f] = fleet;
  } // end fleets loop
  cl_fleets.slot(".Data") = flts;
  return(cl_fleets);
}


