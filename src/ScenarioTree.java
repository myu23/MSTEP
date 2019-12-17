/**
 * This class contains the scenario tree data structure for mutli-stage stochastic programming problem
 * @author Miao Yu
 * @version 1.0
 *
 * Tree structure:
 *      root --------stage 0 (depth 0)
 *     (node0)
 *      /  \
 *   node1 node2 ----stage 1 (depth 1)
 *   .........
 *
 *
 */


import java.util.ArrayList;

public class ScenarioTree {
    public ArrayList<Node> nodeList;// list of nodes
    public ArrayList<ArrayList<Integer>> nodeByDepth; //list of nodes for each depth
    public int maxDepth; //maximum depth (stages) of the tree
    /**
     * Node structure
     */
    public class Node{
        public int index;
        public int depth;
        public int parent;
        public ArrayList<Integer> children;
        public double[] vals; // stochastic data stored in array
        public double probability;

        public Node(int d, int p){
            this.depth = d;
            this.parent = p;
            this.children = new ArrayList<>();
        }

        public void setVals(double[] v){
            this.vals = v;
        }

        public void setProb(double pr){
            this.probability = pr;
        }
    }

    public void addNode(Node n){
        n.index = nodeList.size();
        this.nodeList.add(n);
        if(n.depth != 0){
            this.nodeList.get(n.parent).children.add(n.index);
        }

    }

    public int size(){
        return this.nodeList.size();
    }
    public int nSenario() {return this.nodeByDepth.get(maxDepth).size();}

    /*
    build a tree with (depth+1) levels
     */
    public ScenarioTree(int depth, int nChildPerNode){
        this.nodeList = new ArrayList<>();
        this.nodeByDepth = new ArrayList<>(depth);
        this.maxDepth = depth;

        Node root = new Node(0, -1);
        root.setProb(1.0);
        addNode(root);
        nodeByDepth.add(new ArrayList<>());
        nodeByDepth.get(0).add(root.index);

        ArrayList<Integer> prev = nodeByDepth.get(0);
        for(int d = 0; d < depth; d++){
            nodeByDepth.add(new ArrayList<>());
            for (int i = 0; i < nodeByDepth.get(d).size(); i++){
                Node p = nodeList.get(nodeByDepth.get(d).get(i));
                for(int j = 0; j < nChildPerNode; j++){
                    Node child = new Node(p.depth + 1, p.index);
                    child.setProb(Math.pow(1.0/nChildPerNode, child.depth));
                    addNode(child);
                    nodeByDepth.get(d+1).add(child.index);
                }
            }
        }
    }

    public static void test(){
        ScenarioTree st = new ScenarioTree(10, 3);
        System.out.println(st.nodeByDepth.get(0));
        System.out.println(st.nodeByDepth.get(1));
        System.out.println(st.nodeByDepth.get(10));
        System.out.println(st.nodeList.get(29524).children);
    }
    public static void main(String[] args){
        test();
    }

}
