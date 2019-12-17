public class testDrive {

    public static void main(String[] args){
        Data d = new Data();
        d.loadData(6, 15);
        d.loadSamples(1);
        d.buildTree();

//        DEmodel de = new DEmodel(d, d.tree);
//        de.main();
        SDDiP sddip = new SDDiP(d);
        sddip.solve();
    }

}
