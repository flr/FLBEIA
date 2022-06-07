#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
S4 fill_flcatches(S4 fl,
                  S4 cobj,
                  CharacterVector st,
                  int mt_idx) {

  //the object we want to feel is fl (a fleet)


    S4 x = clone(fl); // a fleet

     List met = x.slot("metiers");             //  metiers

     S4 mt = met[mt_idx-1];

       S4 catches = mt.slot("catches");
       CharacterVector catchNames = catches.slot("names");
       List CL = catches.slot(".Data");   // the FLCatch Data as a list

      if(is_true(any(in(st, catchNames)))) {

         IntegerVector pos = match(st, catchNames) -1;     // this is an R function, so references from 1
         int pos2 = pos[0];      // trick to get back to int from IntegerVector

        S4 C = CL[pos2];                               // The specific stock data

        S4 flqn = C.slot("landings.n");                // just used to get the right dimensions
        NumericVector dimn = flqn.attr("dim");          // save the dimensions attribute

        S4 flq = C.slot("landings");                // just used to get the right dimensions
        NumericVector dim = flq.attr("dim");          // save the dimensions attribute


    //the object from where to fill fl is cobj
    S4 y = cobj;

    S4 yLn = y.slot("landings.n");     // landings.n
    NumericVector yLn_dat = yLn.slot(".Data");

    S4 yDn = y.slot("discards.n");      // discards.n
    NumericVector yDn_dat = yDn.slot(".Data");

    S4 yL = y.slot("landings");     // landings.n
    NumericVector yL_dat = yL.slot(".Data");

    S4 yD = y.slot("discards");      // discards.n
    NumericVector yD_dat = yD.slot(".Data");

    S4 yLwt = y.slot("landings.wt");     // landings.wt
    NumericVector yLwt_dat = yLwt.slot(".Data");

    S4 yDwt = y.slot("discards.wt");      // discards.wt
    NumericVector yDwt_dat = yDwt.slot(".Data");

    S4 yLsel = y.slot("landings.sel");     // landings.sel
    NumericVector yLsel_dat = yLsel.slot(".Data");

    S4 yDsel = y.slot("discards.sel");      // discards.sel
    NumericVector yDsel_dat = yDsel.slot(".Data");

    S4 yq = y.slot("catch.q");     // landings.sel
    NumericVector yq_dat = yq.slot(".Data");

    S4 yalpha = y.slot("alpha");     // landings.sel
    NumericVector yalpha_dat = yalpha.slot(".Data");

    S4 ybeta = y.slot("beta");     // landings.sel
    NumericVector ybeta_dat = ybeta.slot(".Data");

    S4 yP = y.slot("price");     // price
    NumericVector yP_dat = yP.slot(".Data");



    // Return the dimensions to the objects
    yLn_dat.attr("dim") = dimn;
    yDn_dat.attr("dim") = dimn;
    yL_dat.attr("dim") = dim;
    yD_dat.attr("dim") = dim;
    yLwt_dat.attr("dim") = dimn;
    yDwt_dat.attr("dim") = dimn;
    yLsel_dat.attr("dim") = dimn;
    yDsel_dat.attr("dim") = dimn;
    yq_dat.attr("dim") = dimn;
    yalpha_dat.attr("dim") = dimn;
    ybeta_dat.attr("dim") = dimn;
    yP_dat.attr("dim") = dimn;

   //filling needed slots

    S4 Ln_slot = C.slot("landings.n");
    Ln_slot.slot(".Data") = yLn_dat;
    C.slot("landings.n") = Ln_slot;

    S4 L_slot = C.slot("landings");
    L_slot.slot(".Data") = yL_dat;
    C.slot("landings") = L_slot;

    S4 Dn_slot = C.slot("discards.n");
    Dn_slot.slot(".Data") = yDn_dat;
    C.slot("discards.n") = Dn_slot;

    S4 D_slot = C.slot("discards");
    D_slot.slot(".Data") = yD_dat;
    C.slot("discards") = D_slot;

    S4 Lwt_slot = C.slot("landings.wt");
    Lwt_slot.slot(".Data") = yLwt_dat;
    C.slot("landings.wt") = Lwt_slot;

    S4 Dwt_slot = C.slot("discards.wt");
    Dwt_slot.slot(".Data") = yDwt_dat;
    C.slot("discards.wt") = Dwt_slot;

    S4 Lsel_slot = C.slot("landings.sel");
    Lsel_slot.slot(".Data") = yLsel_dat;
    C.slot("landings.sel") = Lsel_slot;

    S4 Dsel_slot = C.slot("discards.sel");
    Dsel_slot.slot(".Data") = yDsel_dat;
    C.slot("discards.sel") = Dsel_slot;

    S4 q_slot = C.slot("catch.q");
    q_slot.slot(".Data") = yq_dat;
    C.slot("catch.q") = q_slot;

    S4 alpha_slot = C.slot("alpha");
    alpha_slot.slot(".Data") = yalpha_dat;
    C.slot("alpha") = alpha_slot;

    S4 beta_slot = C.slot("beta");
    beta_slot.slot(".Data") = ybeta_dat;
    C.slot("beta") = beta_slot;

    S4 P_slot = C.slot("price");
    P_slot.slot(".Data") = yP_dat;
    C.slot("price") = P_slot;


        // Return to FLcatches //

        CL[pos2] = C; // Return to CL
        catches.slot(".Data") = CL; // Return CL to catches

       }
      mt.slot("catches") = catches; // return FLCatches to the metier
      met[mt_idx-1] = mt;                  // return metier to the list of metier

    x.slot("metiers") = met; // return metier to fleet

    //fl = x;   // return fleet to the fleet object

  return(x);


}
