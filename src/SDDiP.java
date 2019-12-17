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
        this.iterCounter = 0;
        this.maxIter = 1;
        this.M = 40;
        this.initPsi = new double[data.T+1];
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

        Psi[][] psi = new Psi[maxIter+1][data.T+1];
        for(int i = 0; i < maxIter+1; i++){
            for(int t = 0; t < data.T+1; t++){
                psi[i][t] = new Psi(initPsi[t]);
            }
        }

        Xi[] xi = new Xi[M];

        VarV[][][] varV = new VarV[maxIter][data.scenarioPerStage][data.T+1];
        VarPi[][][] varPi = new VarPi[maxIter][data.scenarioPerStage][data.T+1];


        while(iterCounter < maxIter){
            // Sampling steps
            for(int k = 0; k < M; k++){
                xi[k] = new Xi();
            }

            // forward steps
            double[] u = new double[M];
            for(int k = 0; k < M; k++){
                for(int t = 0; t < data.T+1; t++){
                    ForwardProblem fp = new ForwardProblem();
                    if(t == 0){
                        fp.solveRoot(psi[iterCounter][0]);
                    }else{
                        fp.solve(varX[iterCounter][k][t-1], psi[iterCounter][t], xi[k], t);
                    }
                    varX[iterCounter][k][t] = fp.x;
                    varY[iterCounter][k][t] = fp.y;
                    varZ[iterCounter][k][t] = fp.z;
                    varTheta[iterCounter][k][t] = fp.theta;
                    u[k] += fp.obj/Math.pow(1.0/data.scenarioPerStage, t);
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

//            used to debug
//            if(iterCounter+1 == maxIter)
//                return;
            // backward steps
            for(int t = data.T; t >= 1; t--){
                for(int k = 0; k < M; k++) {
                    double[] v = new double[data.scenarioPerStage];
                    double[][][] pi = new double[data.scenarioPerStage][data.nBranch][data.maxAddedLine];
                    for (int j = 0; j < data.scenarioPerStage; j++) {
                        BackwardProblem bp = new BackwardProblem();
                        bp.solveBender(varX[iterCounter][k][t - 1], psi[iterCounter+1][t], xi[k], t);
                        varV[iterCounter][j][t] = bp.v;
                        v[j] = bp.v.v;
                        varPi[iterCounter][j][t] = bp.pi;
                        pi[j] = bp.pi.pi;
                    }
                    psi[iterCounter+1][t-1].vList.add(v);
                    psi[iterCounter+1][t-1].piList.add(pi);
                }
            }

            // update lowerbound
            ForwardProblem fp = new ForwardProblem();
            fp.solveRoot(psi[iterCounter+1][0]);
            System.out.println(psi[iterCounter+1][0].vList.size());
            System.out.println(psi[iterCounter+1][0].piList.size());

            lowerbound = fp.obj;

            iterCounter++;
            System.out.println("SDDiP-Iter:"+iterCounter+", lowerbound:"+lowerbound+", upperbound:"+upperbound);
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
            nodes[data.T] = data.tree.nodeByDepth.get(data.T).get(0)+rand;
            for(int t = data.T-1; t >= 0; t--){
                nodes[t] = data.tree.nodeList.get(nodes[t+1]).parent;
            }
        }
    }

    /**
     * Classes for general notations X, Y, Z, Theta, V, Pi ...
     */
    public class VarX{
        double[][] alpha;
        public VarX(double[][] x){
            this.alpha = x;
        }
    }

    public class VarY{
        double[][] flow;
        double[] angle;
        double[] generation;
        double[] demandCurtailment;

        public VarY(double[][] f, double[] a, double[] g, double[] dc){
            this.flow = f;
            this.angle = a;
            this.generation = g;
            this.demandCurtailment = dc;
        }
    }

    public class VarZ{
        double[][] z;

        public VarZ(double[][] z){
            this.z = z;
        }
    }

    public class VarTheta{
        double theta;
        public VarTheta(double t){
            this.theta = t;
        }
    }

    public class Psi{
        public double initBound;
        public ArrayList<double[]> vList;
        public ArrayList<double[][][]> piList;

        public Psi(double bound){
            this.initBound = bound;
            this.vList = new ArrayList<>();
            this.piList = new ArrayList<>();
        }

    }

    public class VarV{
        double v;
        public VarV(double obj){ this.v = obj;}
    }

    public class VarPi{
        double[][] pi;
        public VarPi(double [][] pi){
            this.pi = pi;
        }
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

        public ForwardProblem(){

        }


        public void solve(VarX old_x , Psi psi, Xi xi, int t){
            try {
                GRBEnv env = new GRBEnv();
                GRBModel model = new GRBModel(env);

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
                                expr.addTerm(-psi.piList.get(l)[m][br][k]/((double) data.scenarioPerStage), alpha[br][k]);
                                rhs += psi.vList.get(l)[m]/((double) data.scenarioPerStage);
                            }
                            expr.addTerm(1, theta);
                            model.addConstr(expr, GRB.GREATER_EQUAL, rhs, "StateLink_"+br+"_"+k);

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

                model.optimize();

                //System.out.println("Theta: "+theta.get(GRB.DoubleAttr.X));

                //get solution
                double[][] xsol = model.get(GRB.DoubleAttr.X, alpha);
                double[][] flowsol = model.get(GRB.DoubleAttr.X, flow);
                double[] gsol = model.get(GRB.DoubleAttr.X, generation);
                double[][] zsol = model.get(GRB.DoubleAttr.X, z);
                double[] angsol = model.get(GRB.DoubleAttr.X, angle);
                double[] dcsol = model.get(GRB.DoubleAttr.X, DC);

                this.x = new VarX(xsol);
                this.y = new VarY(flowsol, angsol, gsol, dcsol);
                this.z = new VarZ(zsol);
                this.obj = model.get(GRB.DoubleAttr.ObjVal);
                this.theta = new VarTheta(obj);

                model.dispose();
                env.dispose();



            }catch (GRBException e) {
                System.out.println("Error code: " + e.getErrorCode() + ". " +
                        e.getMessage());
            }


        }


        public void solveRoot(Psi psi){
            try {
                GRBEnv env = new GRBEnv();
                GRBModel model = new GRBModel(env);

                //state variables
                GRBVar[][] alpha = new GRBVar[data.nBranch][data.maxAddedLine];
                for (int br = 0; br < data.nBranch; br++) {
                    for (int k = 0; k < data.maxAddedLine; k++) {
                        String st = "Alpha_" + String.valueOf(br) + "_" + String.valueOf(k);
                        if(k < data.n0[br])
                            alpha[br][k] = model.addVar(1.0, 1.0, 0.0, GRB.BINARY, st);
                        else
                            alpha[br][k] = model.addVar(0, 1.0, 0.0, GRB.BINARY, st);
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
                                expr.addTerm(-psi.piList.get(l)[m][br][k]/((double) data.scenarioPerStage), alpha[br][k]);
                                rhs += psi.vList.get(l)[m]/((double) data.scenarioPerStage);
                            }
                            expr.addTerm(1, theta);
                            model.addConstr(expr, GRB.GREATER_EQUAL, rhs, "StateLink_"+br+"_"+k);

                        }
                    }
                }

                double cost = 0;
                for(int br = 0; br < data.nBranch; br++){
                    for(int k = 0; k < data.n0[br];k++){
                        cost += data.cost[br];
                    }
                }
                GRBVar dummyCost = model.addVar(cost,cost,0,GRB.CONTINUOUS,"dummy_cost");


                //objective
                expr = new GRBLinExpr();
                for(int k = 0; k < data.maxAddedLine; k++){
                    for(int br = 0; br < data.nBranch; br++){
                        expr.addTerm(data.cost[br], alpha[br][k]);
                    }
                }
                expr.addTerm(-1, dummyCost);
                expr.addTerm(1, theta);

                model.setObjective(expr, GRB.MINIMIZE);
                model.write("Root.lp");
                model.optimize();


                //get solution
                double[][] xsol = model.get(GRB.DoubleAttr.X, alpha);

                this.x = new VarX(xsol);
                this.y = null;
                this.z = null;
                this.obj = model.get(GRB.DoubleAttr.ObjVal);
                this.theta = new VarTheta(obj);

                model.dispose();
                env.dispose();



            }catch (GRBException e) {
                System.out.println("Error code: " + e.getErrorCode() + ". " +
                        e.getMessage());
            }


        }
    }

    /**
     * class for backword problem
     */
    public class BackwardProblem{

        public VarV v;
        public VarPi pi;

        public double obj;

        public BackwardProblem(){

        }


        public void solveBender(VarX old_x , Psi psi, Xi xi, int t){
            try {
                GRBEnv env = new GRBEnv();
                GRBModel model = new GRBModel(env);

                //state variables
                GRBVar[][] alpha = new GRBVar[data.nBranch][data.maxAddedLine];
                for (int br = 0; br < data.nBranch; br++) {
                    for (int k = 0; k < data.maxAddedLine; k++) {
                        String st = "Alpha_" + String.valueOf(br) + "_" + String.valueOf(k);
                        if(k < data.n0[br])
                            alpha[br][k] = model.addVar(1.0, 1.0, 0.0, GRB.CONTINUOUS, st);
                        else
                            alpha[br][k] = model.addVar(old_x.alpha[br][k], 1.0, 0.0, GRB.CONTINUOUS, st);
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
                        z[br][k] = model.addVar(0, 1, 0, GRB.CONTINUOUS, "Z-"+br+"-"+k);
                    }
                }

                GRBConstr[][] constrZ = new GRBConstr[data.nBranch][data.maxAddedLine];
                for (int br = 0; br < data.nBranch; br++) {
                    for (int k = 0; k < data.maxAddedLine; k++) {
                        constrZ[br][k] = model.addConstr(old_x.alpha[br][k], GRB.EQUAL, z[br][k],"constrZ-"+br+"-"+k);
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

                model.optimize();



                //get solution
                double[][] pisol = model.get(GRB.DoubleAttr.Pi, constrZ);

                double updatedV = model.get(GRB.DoubleAttr.ObjVal);

                for(int b = 0; b < data.nBranch; b++){
                    for(int k = 0; k < data.maxAddedLine;k++){
                        updatedV -= pisol[b][k]*old_x.alpha[b][k];
                    }
                }
                this.pi = new VarPi(pisol);
                this.v = new VarV(updatedV);

                model.dispose();
                env.dispose();



            }catch (GRBException e) {
                System.out.println("Error code: " + e.getErrorCode() + ". " +
                        e.getMessage());
            }


        }

    }

}
