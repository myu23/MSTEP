import java.util.ArrayList;

public class ScenarioTree {
    public ArrayList<Node> nodeList;

    public class Node{
        public int index;
        public int depth;
        public int parent;
        public ArrayList<Integer> children;
        public double[] vals;

        public Node(int i, int d, int p){
            this.index = i;
            this.depth = d;
            this.parent = p;
            this.children = new ArrayList<>();
        }

        public void setVals(double[] v){
            this.vals = v;
        }
    }

    public void addNode(Node n){
        this.nodeList.add(n);
        if(n.depth != 0){
            this.nodeList.get(n.parent).children.add(n.index);
        }
    }

    public ScenarioTree(int depth, int nChildPerNode){
        this.nodeList = new ArrayList<>();
        Node root = new Node(0, 0, -1);
        nodeList.add(root);
        ArrayList<Integer> prev = new ArrayList<>();
        prev.add(0);
        while(prev.size() != 0){
            if(nodeList.get(prev.get(0)).depth == depth)
                break;
            ArrayList<Integer> curr = new ArrayList<>();
            for (int i = 0; i < prev.size(); i++){
                Node p = nodeList.get(prev.get(i));

            }
        }


    }

}
