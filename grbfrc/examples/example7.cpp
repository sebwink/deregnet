#include <gurobi_c++.h>
#include <grbfrc.h>

#include <cmath>
#include <fstream>

using FMILP = grbfrc::FMILP;

class myCallback: public GRBCallback       
// adapted from Gurobi example callback_c++.cpp
{
  public:

    double lastiter;
    double lastnode;
    int numvars;
    GRBVar* vars;
    std::ofstream* logfile;

    myCallback(int xnumvars, GRBVar* xvars, std::ofstream* xlogfile)
     {
      lastiter = lastnode = -GRB_INFINITY;
      numvars = xnumvars;
      vars = xvars;
      logfile = xlogfile;
     }

  protected:

    void callback () 
     {
      try 
       {
        if (where == GRB_CB_PRESOLVE) 
         {
          // Presolve callback
          int cdels = getIntInfo(GRB_CB_PRE_COLDEL);
          int rdels = getIntInfo(GRB_CB_PRE_ROWDEL);
          if (cdels || rdels) 
           {
            std::cout << cdels << " columns and " << rdels
                      << " rows are removed" << std::endl;
           }
         }
        else if (where == GRB_CB_SIMPLEX) 
         {
          // Simplex callback
          double itcnt = getDoubleInfo(GRB_CB_SPX_ITRCNT);
          if (itcnt - lastiter >= 100)
           {
            lastiter = itcnt;
            double obj = getDoubleInfo(GRB_CB_SPX_OBJVAL);
            int ispert = getIntInfo(GRB_CB_SPX_ISPERT);
            double pinf = getDoubleInfo(GRB_CB_SPX_PRIMINF);
            double dinf = getDoubleInfo(GRB_CB_SPX_DUALINF);
            char ch;
            if (ispert == 0)      ch = ' ';
            else if (ispert == 1) ch = 'S';
            else                  ch = 'P';
            std::cout << itcnt << " " << obj << ch << " "
                      << pinf << " " << dinf << std::endl;
           }
         }
        else if (where == GRB_CB_MIP)
         {
          // General MIP callback
          double nodecnt = getDoubleInfo(GRB_CB_MIP_NODCNT);
          double objbst = getDoubleInfo(GRB_CB_MIP_OBJBST);
          double objbnd = getDoubleInfo(GRB_CB_MIP_OBJBND);
          int solcnt = getIntInfo(GRB_CB_MIP_SOLCNT);
          if (nodecnt - lastnode >= 100)
           {
            lastnode = nodecnt;
            int actnodes = (int) getDoubleInfo(GRB_CB_MIP_NODLFT);
            int itcnt = (int) getDoubleInfo(GRB_CB_MIP_ITRCNT);
            int cutcnt = getIntInfo(GRB_CB_MIP_CUTCNT);
            std::cout << nodecnt << " " << actnodes << " " << itcnt
                      << " " << objbst << " " << objbnd << " "
                      << solcnt << " " << cutcnt << std::endl;
           }
         }
        else if (where == GRB_CB_MIPSOL)
         {
          // MIP solution callback
          int nodecnt = (int) getDoubleInfo(GRB_CB_MIPSOL_NODCNT);
          double obj = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
          int solcnt = getIntInfo(GRB_CB_MIPSOL_SOLCNT);
          double* x = getSolution(vars, numvars);
          std::cout << "**** New solution at node " << nodecnt
                    << ", obj " << obj << ", sol " << solcnt
                    << ", x[0] = " << x[0] << " ****" << std::endl;
          delete[] x;
         }
        else if (where == GRB_CB_MIPNODE)
         {
          // MIP node callback
          std::cout << "**** New node ****" << std::endl;
          if (getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
           {
            double* x = getNodeRel(vars, numvars);
            setSolution(vars, x, numvars);
            delete[] x;
           }
         } 
        else if (where == GRB_CB_MESSAGE) 
         {
          // Message callback
          std::string msg = getStringInfo(GRB_CB_MSG_STRING);
          *logfile << msg;
         }
       }
      catch (GRBException e)
       {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
       }
      catch (...) 
       {
        std::cout << "Error during callback" << std::endl;
       }
     }
  };


int main()
 {
  // define model as usual
  GRBEnv env;
  FMILP model(env);
  GRBVar x = model.addVar(0.0, GRB_INFINITY, GRB_CONTINUOUS, "x");
  GRBVar y = model.addVar(0.0, GRB_INFINITY, GRB_CONTINUOUS, "y");
  GRBVar z = model.addVar(0.0, GRB_INFINITY, GRB_CONTINUOUS, "z");
  GRBVar u = model.addVar(0.0, 1.0, GRB_BINARY, "u");
  GRBVar v = model.addVar(0.0, 1.0, GRB_BINARY, "v");
  model.update();
  model.addConstr( x + 3*y + z <= 10 );
  model.addConstr( 3*x + 2*y <= 7 );
  model.addConstr( u + v >= 1 );
  GRBLinExpr objNum { 2*x + y + 3*z + u + v + 6 };
  model.setObjNumerator(objNum);
  GRBLinExpr objDenom { x + 3*z + v + 4 };
  model.setObjDenominator(objDenom);
  model.setObjSense(GRB_MAXIMIZE);
  // define callback
  GRBVar* vars = model.getVars();
  std::ofstream logfile("cb.log");
  myCallback cb { 5, vars, &logfile };
  // run Dinkelbach's algorithm with callback
  model.runDinkelbach(cb);
  model.printSolution();
 }
