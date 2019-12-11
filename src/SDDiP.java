/**
 * This class contains the SDDiP framework that can be extended to solve other problems.
 * Implement Algorithm 2 in Stochatic Dual Dynamic Integer Programming http://www.optimization-online.org/DB_FILE/2016/05/5436.pdf
 * Stopping criteria is set as maximum iteration reached
 *
 * @author Miao
 * @version 0.1
 */


import gurobi.*;

import java.util.ArrayList;
import java.util.Random;


public class SDDiP {
    private double lowerbound;
    private double upperbound;

    private int iterCounter; //iteration counter
    private int maxIter; //stop criteria by iteration
    private int M;// sample M scenarios (no less than 30)



    private double[] initPsi; //initial under-approximation with dimension of nStage
    private Data data;
    private Random random;





    /**
     * class constructor with input data
     * @param data
     */
    public SDDiP(Data data){
        this.data = data;
        this.lowerbound = Double.MIN_VALUE;
        this.upperbound = Double.MAX_VALUE;
        this.iterCounter = 1;
        this.maxIter = 50;
        this.M = 40;
        this.initPsi = new double[data.T];
        this.random = new Random(1);
    }

    /**
     * manually set initial under-approxiamtion
     * @param psi
     */
    public void setInitPsi(double[] psi){
        this.initPsi = psi;
    }

    /**
     * manually set stop criteria
     * @param niter
     */
    public void setMaxIter(int niter){
        this.maxIter = niter;
    }

    /**
     * implement sddip framework
     */
    public void solve(){

        VarX[][][] varX = new VarX[maxIter][M][data.T+1];
        VarY[][][] varY = new VarY[maxIter][M][data.T+1];
        VarZ[][][] varZ = new VarZ[maxIter][M][data.T+1];
        VarTheta[][][] varTheta = new VarTheta[maxIter][M][data.T+1];

        Psi[][] psi = new Psi[M][data.T+1];
        Xi[] xi = new Xi[M];

        VarV[][][] varV = new VarV[maxIter][M][data.T+1];
        VarPi[][][] varPi = new VarPi[maxIter][M][data.T+1];


        while(iterCounter < maxIter){
            // Sampling steps
            for(int k = 0; k < M; k++){
                xi[k] = new Xi();
            }

            // forward steps
            double[] u = new double[M];
            for(int k = 0; k < M; k++){
                for(int t = 1; t < data.T+1; t++){
                    ForwardProblem fp = new ForwardProblem(varX[iterCounter][k][t-1], psi[iterCounter][t], xi[k], t);
                    fp.solve();
                    varX[iterCounter][k][t] = fp.x;
                    varY[iterCounter][k][t] = fp.y;
                    varZ[iterCounter][k][t] = fp.z;
                    varTheta[iterCounter][k][t] = fp.theta;
                    u[k] += fp.obj;
                }
            }

            // Statistical upper bound update
            double muh = 0;
            for(double d : u) muh += d;
            muh = muh/M;

            double sigma2 = 0;
            for(double d : u) sigma2 += (d - muh)*(d - muh);
            sigma2 = sigma2 / (M-1);

            upperbound = muh + 1.96 * Math.pow(sigma2/M, 0.5);


            // backward steps
            for(int t = data.T+1; t > 1; t--){
                for(int k = 0; k < M; k++){
                    for(int j = 0; j < data.scenarioPerStage; j++){
                        BackwardProblem bp = new BackwardProblem(varX[iterCounter][k][t-1], psi[iterCounter+1][t], xi[k], t);
                        bp.solve();
                        varV[iterCounter][k][t] = bp.v;
                        varPi[iterCounter][k][t] = bp.pi;

                    }
                }
            }

        }
    }

    /**
     * class used to store stochastic data
     */
    public class Xi{
        public int[] nodes;


        public Xi(){
            int nSenario = data.tree.nSenario();
            int rand = random.nextInt(nSenario);
            nodes = new int[data.T+1];
            nodes[data.T] = rand;
            for(int t = data.T-1; t >= 0; t--){
                nodes[t] = data.tree.nodeList.get(nodes[t+1]).parent;
            }
        }
    }

    public class VarX{
        double[][] alpha;
    }

    public class VarY{
        double[][] flow;
        double[] angle;
        double[] generation;
        double[] demandCurtailment;
    }

    public class VarZ{
        double[] z;
    }

    public class VarTheta{

    }

    public class Psi{
        public double initBound;
        public ArrayList<Double[]> vList;
        public ArrayList<Double[][][]> piList;

    }

    public class VarV{

    }

    public class VarPi{

    }
    /**
     * class for forward problem
     */
    public class ForwardProblem{
        public VarX x;
        public VarY y;
        public VarZ z;
        public VarTheta theta;
        public double obj;
        public GRBModel model;

        public ForwardProblem(VarX old_x , Psi psi, Xi xi, int t){
            try {
                GRBEnv env = new GRBEnv();
                model = new GRBModel(env);

                //state variables
                GRBVar[][] alpha = new GRBVar[data.nBranch][data.maxAddedLine];
                for (int br = 0; br < data.nBranch; br++) {
                    for (int k = 0; k < data.maxAddedLine; k++) {
                        String st = "Alpha_" + String.valueOf(br) + "_" + String.valueOf(k);
                        if(k < data.n0[br])
                            alpha[br][k] = model.addVar(1.0, 1.0, 0.0, GRB.BINARY, st);
                        else
                            alpha[br][k] = model.addVar(old_x.alpha[br][k], 1.0, 0.0, GRB.BINARY, st);
                    }
                }
                //obj var
                GRBVar theta = model.addVar(psi.initBound, Double.POSITIVE_INFINITY, 0, GRB.CONTINUOUS, "theta");

                GRBLinExpr expr;
                //state constraint
                for (int br = 0; br < data.nBranch; br++) {
                    for (int k = 1; k < data.maxAddedLine; k++) {
                        expr = new GRBLinExpr();
                        expr.addTerm(1.0, alpha[br][k]);
                        expr.addTerm(-1.0, alpha[br][k-1]);
                        model.addConstr(expr, GRB.LESS_EQUAL, 0, "StateLink_"+br+"_"+k);
                    }
                }

                // add generated cuts
                for(int l = 0; l < psi.vList.size(); l++){
                    expr = new GRBLinExpr();
                    for (int br = 0; br < data.nBranch; br++) {
                        for (int k = 1; k < data.maxAddedLine; k++) {
                            expr = new GRBLinExpr();
                            double rhs = 0;
                            for(int m = 0; m < data.scenarioPerStage; m++) {
                                expr.addTerm(-psi.piList.get(l)[m][br][k], alpha[br][k]);
                                rhs += psi.vList.get(l)[m];
                            }
                            expr.addTerm(1, theta);
                            model.addConstr(expr, GRB.LESS_EQUAL, rhs, "StateLink_"+br+"_"+k);

                        }
                    }
                }

                // Create copys of local variables
                GRBVar[] generation = new GRBVar[data.nBus];
                for(int b = 0; b < data.nBus; b++)
                    generation[b] = model.addVar(0, data.gen_level[b], 0, GRB.CONTINUOUS, "G_"+b);



                GRBVar[][] flow = new GRBVar[data.nBranch][data.maxAddedLine];
                for(int b = 0; b < data.nBranch; b++)
                    for(int k = 0; k < data.maxAddedLine; k++)
                        flow[b][k] = model.addVar(-data.barF[b], data.barF[b], 0, GRB.CONTINUOUS, "Flow_"+b+"_"+k);

                GRBVar[] angle = new GRBVar[data.nBus];
                for(int b = 0; b < data.nBus; b++)
                    angle[b] = model.addVar(-Math.PI/2, Math.PI/2, 0, GRB.CONTINUOUS, "Theta_"+b);

                GRBVar[] DC = new GRBVar[data.nBus];
                for(int b = 0; b < data.nBus; b++) {
                    DC[b] = model.addVar(0, data.tree.nodeList.get(xi.nodes[t]).vals[b], 0, GRB.CONTINUOUS, "DC_" + b);
                }

                //add local constraints

                // connectivity
                // Sf + g + DC = Demand(s)

                for(int b = 0; b < data.nBus; b++){
                    expr = new GRBLinExpr();
                    for(int br = 0; br < data.nBranch; br++){
                        if(data.to[br] == b){
                            for(int k = 0; k < data.maxAddedLine; k++)
                                expr.addTerm(1.0, flow[br][k]);
                        }else if(data.from[br] == b){
                            for(int k = 0; k < data.maxAddedLine; k++)
                                expr.addTerm(-1.0, flow[br][k]);
                        }
                    }
                    expr.addTerm(1.0, generation[b]);
                    expr.addTerm(1.0, DC[b]);
                    model.addConstr(expr,GRB.EQUAL,data.tree.nodeList.get(xi.nodes[t]).vals[b], "FlowBalance_"+b);
                }

                //initial valid
                //sum_flow - gamma_{ij}*n0_{ij}*(theta_i - theta_j) = 0
                for(int br = 0; br < data.nBranch; br++){
                    expr = new GRBLinExpr();
                    for(int k = 0; k < data.n0[br]; k++){
                        expr.addTerm(1.0, flow[br][k]);
                    }
                    expr.addTerm(-100.0 / data.reactance[br] * data.n0[br], angle[data.from[br]]);
                    expr.addTerm(+100.0 / data.reactance[br] * data.n0[br], angle[data.to[br]]);
                    model.addConstr(expr, GRB.EQUAL, 0, "FlowValidInit_"+br);
                }

                //new valid
                //flow - gamma_{ij}(theta_i-theta_j) <= M(1-alpha)
                //flow - gamma_{ij}(theta_i-theta_j) >= -M(1-alpha)
                for(int br = 0; br < data.nBranch; br++) {
                    for (int k = data.n0[br]; k < data.maxAddedLine; k++) {
                        expr = new GRBLinExpr();
                        expr.addTerm(1.0, flow[br][k]);
                        expr.addTerm(-100.0 / data.reactance[br], angle[data.from[br]]);
                        expr.addTerm(+100.0 / data.reactance[br], angle[data.to[br]]);
                        model.addConstr(expr, GRB.LESS_EQUAL, data.M * (1-old_x.alpha[br][k]), "FlowValidNew_" + br);

                        expr = new GRBLinExpr();
                        expr.addTerm(1.0, flow[br][k]);
                        expr.addTerm(-100.0 / data.reactance[br], angle[data.from[br]]);
                        expr.addTerm(+100.0 / data.reactance[br], angle[data.to[br]]);
                        model.addConstr(expr, GRB.GREATER_EQUAL, -data.M*(1-old_x.alpha[br][k]), "FlowValidNew_" + br);
                    }
                }

                //flow capacity
                for(int br = 0; br < data.nBranch; br++) {
                    for (int k = 0; k < data.maxAddedLine; k++) {
                        flow[br][k].set(GRB.DoubleAttr.UB, data.barF[br]*old_x.alpha[br][k]);
                        flow[br][k].set(GRB.DoubleAttr.LB, -data.barF[br]*old_x.alpha[br][k]);
                    }
                }



                GRBVar[][] z = new GRBVar[data.nBranch][data.maxAddedLine];
                for (int br = 0; br < data.nBranch; br++) {
                    for (int k = 0; k < data.maxAddedLine; k++) {
                        z[br][k] = model.addVar(old_x.alpha[br][k], old_x.alpha[br][k], 0, GRB.CONTINUOUS, "Z-"+br+"-"+k);
                    }
                }



                //objective

                expr = new GRBLinExpr();
                for(int k = 0; k < data.maxAddedLine; k++){
                    for(int br = 0; br < data.nBranch; br++){
                        expr.addTerm(Math.pow(data.discountRate, t) * data.cost[br], alpha[br][k]);
                        expr.addTerm(-Math.pow(data.discountRate, t) * data.cost[br], z[br][k]);
                    }
                }
                for(int b = 0; b < data.nBus; b++){
                    expr.addTerm(Math.pow(data.discountRate, t) * data.dcPenalty, DC[b]);
                }
                expr.addTerm(1, theta);

                model.setObjective(expr, GRB.MINIMIZE);



            }catch (GRBException e) {
                System.out.println("Error code: " + e.getErrorCode() + ". " +
                        e.getMessage());
            }


        }


        public void solve(){
            try{
                model.optimize();

            }catch(GRBException e){
                System.out.println("Error code: " + e.getErrorCode() + ". " +
                        e.getMessage());
            }
        };
    }

    /**
     * class for backword problem
     */
    public class BackwardProblem{

        public VarV v;
        public VarPi pi;

        public double obj;

        public BackwardProblem(VarX x , Psi psi, Xi xi, int t){

        }


        public void solve(){};
    }

}
