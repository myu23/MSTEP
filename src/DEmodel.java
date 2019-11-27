/**
 * This class contains the solver to the deterministic equivalent of the multistage
 * transmission extension planing problem.
 *
 *
 */

import gurobi.*;
import java.io.*;
import java.util.Arrays;

public class DEmodel {
    public Data data;
    public ScenarioTree tree;

    public DEmodel(Data d, ScenarioTree t){
        this.data = d;
        this.tree = d.tree;
    }
    public void main(){
        try {
            GRBEnv env = new GRBEnv();
            GRBModel model = new GRBModel(env);

            // Create 3-D array of state variables

            GRBVar[][][] alpha = new GRBVar[data.nBranch][data.maxAddedLine][tree.size()];

            for (int br = 0; br < data.nBranch; br++) {
                for (int k = 0; k < data.maxAddedLine; k++) {
                    for (int i = 0; i < tree.size(); i++) {
                        String st = "Alpha_" + String.valueOf(br) + "_" + String.valueOf(k)
                                + "_" + String.valueOf(i);
                        if(k < data.n0[br])
                            alpha[br][k][i] = model.addVar(1.0, 1.0, 0.0, GRB.BINARY, st);
                        else
                            alpha[br][k][i] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, st);
                    }
                }
            }

            // Create copys of local variables
            GRBVar[][] generation = new GRBVar[data.nBus][tree.size()];
            for(int b = 0; b < data.nBus; b++)
                for(int i = 0; i < tree.size(); i++)
                    generation[b][i] = model.addVar(0, data.gen_level[b], 0, GRB.CONTINUOUS, "G_"+b+"_"+i);



            GRBVar[][][] flow = new GRBVar[data.nBranch][data.maxAddedLine][tree.size()];
            for(int b = 0; b < data.nBranch; b++)
                for(int k = 0; k < data.maxAddedLine; k++)
                    for(int i = 0; i < tree.size(); i++)
                        flow[b][k][i] = model.addVar(-data.barF[b], data.barF[b], 0, GRB.CONTINUOUS, "Flow_"+b+"_"+k+"_"+i);

            GRBVar[][] theta = new GRBVar[data.nBus][tree.size()];
            for(int b = 0; b < data.nBus; b++)
                for(int i = 0; i < tree.size(); i++)
                    theta[b][i] = model.addVar(-Math.PI/2, Math.PI/2, 0, GRB.CONTINUOUS, "Theta_"+b+"_"+i);

            GRBVar[][] DC = new GRBVar[data.nBus][tree.size()];
            for(int b = 0; b < data.nBus; b++) {
                for (int i = 0; i < tree.size(); i++){
                    DC[b][i] = model.addVar(0, 10 * data.load[b], 0, GRB.CONTINUOUS, "DC_" + b + "_" + i);
                }
            }
                // Add constraints

            GRBLinExpr expr;
            //State var constraints
            for(int s = 1; s < tree.size(); s++) {
                int p = tree.nodeList.get(s).parent;
                for (int br = 0; br < data.nBranch; br++) {
                    for (int k = 0; k < data.maxAddedLine; k++) {
                        expr = new GRBLinExpr();
                        expr.addTerm(1.0, alpha[br][k][s]);
                        expr.addTerm(-1.0, alpha[br][k][p]);
                        model.addConstr(expr, GRB.GREATER_EQUAL, 0, "StateLink_"+br+"_"+k+"_"+s);
                    }
                }
            }

            for(int s = 0; s < tree.size(); s++) {
                for (int br = 0; br < data.nBranch; br++) {
                    for (int k = 1; k < data.maxAddedLine; k++) {
                        expr = new GRBLinExpr();
                        expr.addTerm(1.0, alpha[br][k][s]);
                        expr.addTerm(-1.0, alpha[br][k-1][s]);
                        model.addConstr(expr, GRB.LESS_EQUAL, 0, "StateLink_"+br+"_"+k+"_"+s);
                    }
                }
            }

            //scenario based constraints
            //omit root node as we have no local constraints for root node
            for(int s = 1; s < tree.size(); s++){
                int p = tree.nodeList.get(s).parent;

                // connectivity
                // Sf + g + DC = Demand(s)
                for(int b = 0; b < data.nBus; b++){
                    expr = new GRBLinExpr();
                    for(int br = 0; br < data.nBranch; br++){
                        if(data.to[br] == b){
                            for(int k = 0; k < data.maxAddedLine; k++)
                                expr.addTerm(1.0, flow[br][k][s]);
                        }else if(data.from[br] == b){
                            for(int k = 0; k < data.maxAddedLine; k++)
                                expr.addTerm(-1.0, flow[br][k][s]);
                        }
                    }
                    expr.addTerm(1.0, generation[b][s]);
                    expr.addTerm(1.0, DC[b][s]);
                    model.addConstr(expr,GRB.EQUAL, 2*tree.nodeList.get(s).vals[b], "FlowBalance_"+b+"_"+s);
                }

                //initial valid
                //sum_flow - gamma_{ij}*n0_{ij}*(theta_i - theta_j) = 0
                for(int br = 0; br < data.nBranch; br++){
                    expr = new GRBLinExpr();
                    for(int k = 0; k < data.maxAddedLine; k++){
                        expr.addTerm(1.0, flow[br][k][s]);
                    }
                    expr.addTerm(-100.0 / data.reactance[br] * data.n0[br], theta[data.from[br]][s]);
                    expr.addTerm(+100.0 / data.reactance[br] * data.n0[br], theta[data.to[br]][s]);
                    model.addConstr(expr, GRB.EQUAL, 0, "FlowValidInit_"+br+"_"+s);
                }

                //new valid
                //flow - gamma_{ij}(theta_i-theta_j) <= M(1-alpha)
                //flow - gamma_{ij}(theta_i-theta_j) >= -M(1-alpha)
                for(int br = 0; br < data.nBranch; br++) {
                    for (int k = 0; k < data.maxAddedLine; k++) {
                        expr = new GRBLinExpr();
                        expr.addTerm(1.0, flow[br][k][s]);
                        expr.addTerm(-100.0 / data.reactance[br], theta[data.from[br]][s]);
                        expr.addTerm(+100.0 / data.reactance[br], theta[data.to[br]][s]);
                        expr.addTerm(data.M, alpha[br][k][p]);
                        model.addConstr(expr, GRB.LESS_EQUAL, data.M, "FlowValidNew_" + br + "_" + s);

                        expr = new GRBLinExpr();
                        expr.addTerm(1.0, flow[br][k][s]);
                        expr.addTerm(-100.0 / data.reactance[br], theta[data.from[br]][s]);
                        expr.addTerm(+100.0 / data.reactance[br], theta[data.to[br]][s]);
                        expr.addTerm(-data.M, alpha[br][k][p]);
                        model.addConstr(expr, GRB.GREATER_EQUAL, -data.M, "FlowValidNew_" + br + "_" + s);
                    }
                }

                //flow capacity
                for(int br = 0; br < data.nBranch; br++) {
                    for (int k = 0; k < data.maxAddedLine; k++) {
                        expr = new GRBLinExpr();
                        expr.addTerm(1.0, flow[br][k][s]);
                        expr.addTerm(-data.barF[br], alpha[br][k][p]);
                        model.addConstr(expr, GRB.LESS_EQUAL, 0, "FlowCap_"+br+"_"+k+"_"+s);

                        expr = new GRBLinExpr();
                        expr.addTerm(-1.0, flow[br][k][s]);
                        expr.addTerm(-data.barF[br], alpha[br][k][p]);
                        model.addConstr(expr, GRB.LESS_EQUAL, 0, "FlowCap_"+br+"_"+k+"_"+s);
                    }
                }




            }

            // Objective
            expr = new GRBLinExpr();
            for(int br = 0; br < data.nBranch; br++)
                for(int k = data.n0[br]; k < data.maxAddedLine; k++){
                    expr.addTerm(data.cost[br], alpha[br][k][0]);
                }

            for(int t = 1; t <= data.T; t++){
                for(int s : tree.nodeByDepth.get(t)) {
                    //System.out.println(s +"_"+tree.nodeByDepth.get(t).toString());
                    for(int br = 0; br < data.nBranch; br++)
                        for(int k = data.n0[br]; k < data.maxAddedLine; k++) {
                            expr.addTerm(1.0 / (tree.nodeByDepth.get(t).size())
                                    * Math.pow(data.discountRate, t) * data.cost[br], alpha[br][k][s]);
                            expr.addTerm(-1.0 / (tree.nodeByDepth.get(t).size())
                                            * Math.pow(data.discountRate, t) * data.cost[br],
                                    alpha[br][k][tree.nodeList.get(s).parent]);
                        }

                    for(int b = 0; b < data.nBus; b++){
                        expr.addTerm(1.0 / (tree.nodeByDepth.get(t).size())
                                * Math.pow(data.discountRate, t) * 10000, DC[b][s]);
                    }
                }
            }
            model.setObjective(expr,GRB.MINIMIZE);

            // Optimize model

            model.optimize();

            // Write model to file
            model.write("TEP.lp");

            double[][][] x = model.get(GRB.DoubleAttr.X, alpha);

            System.out.println();
            for (int v = 0; v < tree.size(); v++) {
                System.out.println("treenode "+v);
                for (int i = 0; i < data.nBranch; i++) {
                    for (int j = 0; j < data.maxAddedLine; j++){
                        System.out.print(x[i][j][v]+" ");
                    }
                    System.out.println();
                }
            }

            // Dispose of model and environment
            model.dispose();
            env.dispose();

        } catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " +
                    e.getMessage());
        }
    }

}




