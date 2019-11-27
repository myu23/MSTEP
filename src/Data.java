import java.io.*;
import java.util.Arrays;

public class Data {

    //parameters
    public int T;
    public int scenarioPerStage;
    public double discountRate;
    public int maxAddedLine;

    public int default_T = 10;
    public int default_scenarioPerStage = 5;
    public int M = 1000;
    // branch data
    public int nBranch;
    public int[] n0;
    public int[] from;
    public int[] to;
    public double[] reactance;
    public int[] barF;
    public double[] cost;

    // bus data
    public int nBus;
    public int[] gen_max;
    public int[] gen_level;
    public double[] load;
    public double[][][] stochasticLoad;

    // Scenario Tree
    public ScenarioTree tree;
    /**
     *Default constructor
     */
    public Data(){
        this.T = 5;
        this.scenarioPerStage = 3;
        this.discountRate = 1.0/1.1;
        this.maxAddedLine = 6;
    }
    /**
     *Load data from csv file
     *@param nBus: number of bus
     *@param nBranch: number of branch
     */
    public void loadData(int nBus, int nBranch){
        System.out.println("Start loading topology data...");
        String dir = "./Data/"+nBus+"bus/";
        BufferedReader br = null;
        String line = "";
        String cvsSplitBy = ",";

        //bus
        this.nBus = nBus;
        this.gen_level = new int[nBus];
        this.gen_max = new int[nBus];
        this.load = new double[nBus];
        String busFile = dir+"bus_data.csv";
        try {
            br = new BufferedReader(new FileReader(busFile));
            line = br.readLine();
            for(int i = 0; i < nBus; i++){
                line = br.readLine();
                String[] li = line.split(cvsSplitBy);
                gen_max[i] = Integer.parseInt(li[1]);
                gen_level[i] = Integer.parseInt(li[2]);
                load[i] = Double.parseDouble(li[3]);
            }
            System.out.println("Bus data loaded.");

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        //branch data
        //bus
        this.nBranch = nBranch;
        this.n0 = new int[nBranch];
        this.from = new int[nBranch];
        this.to = new int[nBranch];
        this.reactance = new double[nBranch];
        this.barF = new int[nBranch];
        this.cost = new double[nBranch];
        String branchFile = dir+"branch_data.csv";
        try {
            br = new BufferedReader(new FileReader(branchFile));
            line = br.readLine();
            for(int i = 0; i < nBranch; i++){
                line = br.readLine();
                String[] li = line.split(cvsSplitBy);
                from[i] = Integer.parseInt(li[0])-1;
                to[i] = Integer.parseInt(li[1])-1;
                n0[i] = Integer.parseInt(li[2]);
                reactance[i] = Double.parseDouble(li[3]);
                barF[i] = Integer.parseInt(li[4]);
                cost[i] = Double.parseDouble(li[5]);
            }
            System.out.println("Branch data loaded.");


        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

    }

    public void test(){
        Data d = new Data();
        d.generateInstance(6,15, 1);
        System.out.println(Arrays.toString(d.gen_level));
        System.out.println(Arrays.toString(d.cost));
        System.out.println(Arrays.toString(d.stochasticLoad[1][0]));
        System.out.println(Arrays.toString(d.tree.nodeList.get(4).vals));

    }

    /**
     * load stochastic samples from csv files
     * @param number number of samples to load
     */
    public void loadSamples(int number){
        System.out.println("Start loading samples...");
        String dir = "./Sample/"+nBus+"bus/";
        BufferedReader br = null;
        String line = "";
        String cvsSplitBy = ",";


        String sampleFile = dir+"StochLoad.csv";
        this.stochasticLoad = new double[T][scenarioPerStage][nBus];
        try {
            br = new BufferedReader(new FileReader(sampleFile));
            line = br.readLine();
            // skip first 'number-1' lines
            for(int i = 0; i < number; i++){
                line = br.readLine();
            }
            line = br.readLine();
            String[] li = line.split(cvsSplitBy);
            int count = -1;
            for(int t = 0; t < default_T; t++){
                for(int s = 0; s < default_scenarioPerStage; s++) {
                    for(int b = 0; b < nBus; b++){
                        count++;
                        if(t >= T || s >= scenarioPerStage)
                            continue;
                        else
                            stochasticLoad[t][s][b] = Double.parseDouble(li[count]);
                    }
                }
            }
            System.out.println("Sample loaded.");

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) {
                try {
                    br.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    public void buildTree(){
        this.tree = new ScenarioTree(T, scenarioPerStage);
        for(int t = 0; t < T; t++){
            for(int i = 0; i < tree.nodeByDepth.get(t+1).size(); i++){
                ScenarioTree.Node node = tree.nodeList.get(tree.nodeByDepth.get(t+1).get(i));
                //System.out.println(t+" "+Arrays.toString(stochasticLoad[t][i%scenarioPerStage]));
                node.setVals(stochasticLoad[t][i%scenarioPerStage]);
            }
        }
        System.out.println("Tree built.");
    }

    public void generateInstance(int bus, int branch, int idx){
        loadData(bus, branch);
        loadSamples(idx);
        buildTree();
    }

}
