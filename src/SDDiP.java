/**
 * This class contains the SDDiP framework that can be extended to solve other problems.
 * Implement Algorithm 2 in Stochatic Dual Dynamic Integer Programming http://www.optimization-online.org/DB_FILE/2016/05/5436.pdf
 * Stopping criteria is set as maximum iteration reached
 *
 * @author Miao
 * @version 0.1
 */


import gurobi.*;

import java.util.Random;


public class SDDiP {
    private double lowerbound;
    private double upperbound;

    private int iterCounter; //iteration counter
    private int maxIter; //stop criteria by iteration
    private int M;// sample M scenarios (no less than 30)



    private double[] initPsi; //initial under-approximation with dimension of nStage
    private Data data;





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
    }

    /**
     * manually set initial under-approxiamtion
     * @param psi
     */
    public void setPsi(double[] psi){
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
                        BackwardProblem bp = new BackwardProblem(varX[iterCounter][k][t-1], psi[iterCounter][t], xi[k]);
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
            int total = data.tree.nodeByDepth.get(data.T).size();
            int rand = (int) (Math.random()* total);
            nodes = new int[data.T+1];
            nodes[data.T] = rand;
            for(int t = data.T -1; t >= 0; t--){
                nodes[t] = data.tree.nodeList.get(nodes[t+1]).parent;
            }
        }
    }

    public class VarX{

    }

    public class VarY{

    }

    public class VarZ{

    }

    public class VarTheta{

    }

    public class Psi{

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

        public ForwardProblem(VarX x , Psi psi, Xi xi, int t){

        }


        public void solve(){};
    }

    /**
     * class for backword problem
     */
    public class BackwardProblem{

        public VarV v;
        public VarPi pi;

        public double obj;

        public BackwardProblem(VarX x , Psi psi, Xi xi){

        }


        public void solve(){};
    }

}
