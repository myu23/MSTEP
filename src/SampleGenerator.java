import java.io.*;

public class SampleGenerator {
    public double[][][] stochasticLoad; // dim: T, scenarioPerStage, nBus
    public double[][][] stochasticGen;  // dim: T, scenarioPerStage, nBus
    public Data data;

    /**
     * default constructor
     * @param data
     */
    public SampleGenerator(Data data){
        this.stochasticGen = new double[data.T][data.scenarioPerStage][data.nBus];
        this.stochasticLoad = new double[data.T][data.scenarioPerStage][data.nBus];
        this.data = data;
    }

    public void generate(){
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new File("./Sample/"+data.nBus+"bus/StochLoad.csv"));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        StringBuilder builder = new StringBuilder();

        //column name t-s-bus
        for(int t = 0; t < data.default_T; t++){
            for(int s = 0; s < data.default_scenarioPerStage; s++) {
                for (int b = 0; b < data.nBus; b++) {
                    builder.append(t + "-" + s + "-" + b + ",");
                }
            }
        }
        builder.deleteCharAt(builder.length() - 1);
        builder.append("\n");


        for(int i = 0; i < 1000; i++){

            for(int t = 0; t < data.default_T; t++){
                for(int s = 0; s < data.default_scenarioPerStage; s++) {
                    for(int b = 0; b < data.nBus; b++){
                        stochasticLoad[t][s][b] = data.load[b] * (0.6 + Math.random()*0.4) * Math.pow(1.05,t);
                        builder.append(stochasticLoad[t][s][b]);
                        builder.append(",");
                    }
                }
            }
            builder.deleteCharAt(builder.length() - 1);
            builder.append("\n");
        }
        pw.write(builder.toString());
        pw.close();
    }

    public static void main(String[] args){
        Data data = new Data();
        data.loadData(46,79);
        SampleGenerator sg = new SampleGenerator(data);

        sg.generate();

    }
}





