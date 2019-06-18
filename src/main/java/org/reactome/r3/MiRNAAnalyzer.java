/*
 * Created on Apr 14, 2009
 *
 */
package org.reactome.r3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

public class MiRNAAnalyzer {
    
    public MiRNAAnalyzer() {
    }
    

    /**
     * Load a map from a gene to its miRNAs.
     * @return
     * @throws IOException
     */
    public Map<String, Set<String>> loadGeneToMiRNAs() throws IOException {
        String fileName = "datasets/miRNA/miRNATargets/mirecords_mirtarbase_evidenceFiltered_02JUL2014.txt";
        Map<String, Set<String>> miRNAToTargets = new FileUtility().loadSetMap(fileName);
        return InteractionUtilities.switchKeyValues(miRNAToTargets);
    }
    
    @Test
    public void checkDegrees() throws IOException {
        String dir = "/Users/gwu/Dropbox/Eric_Share/processed_data/miRNA_DB/";
//        String file = dir + "mirecords_mirtarbase_inreactome_070114.txt";
        String file = dir + "mirecords_mirtarbase_evidenceFiltered_02JUL2014.txt";
        FileUtility fu = new FileUtility();
        Map<String, Set<String>> miRNAToTargets = fu.loadSetMap(file);
        for (String miRNA : miRNAToTargets.keySet()) {
            System.out.println(miRNA + "\t" + miRNAToTargets.get(miRNA).size());
        }
    }
    
    @Test
    public void checkTargets() throws Exception {
//    	String dir = "/Users/Gwu/Dropbox/Eric_Share/processed_data/miRNA_Correlation/miR199_miR214/";
//    	String file = dir + "gene_list_miR-199.txt";
    	String dir = "/Users/Gwu/Dropbox/Eric_Share/processed_data/miRNA_Correlation/miR103_miR107/";
    	String file = dir + "gene_list_miR-107.txt";
    	FileUtility fu = new FileUtility();
    	Set<String> miR199Genes = fu.loadInteractions(file);
//    	file = dir + "gene_list_miR-214.txt";
    	file = dir + "gene_list_miR-103a.txt";
    	Set<String> miR214Genes = fu.loadInteractions(file);
    	Set<String> shared = InteractionUtilities.getShared(miR199Genes, miR214Genes);
    	double pvalue = MathUtilities.calculateHypergeometricPValue(11000, // Total genes in the FIN
    																miR199Genes.size(),
    																miR214Genes.size(), 
    																shared.size());
    	System.out.println("miR199 genes: " + miR199Genes.size());
    	System.out.println(miR199Genes);
    	System.out.println("miR214 genes: " + miR214Genes.size());
    	System.out.println(miR214Genes);
    	System.out.println("Shared: " + shared.size());
    	System.out.println(shared);
    	System.out.println("P-value: " + pvalue);
    }
    
    @Test
    public void checkMiRecordsFile() throws IOException {
        String fileName = "datasets/miRNA/miRecords/miRecords_version3.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> fis = new HashSet<String>();
        Set<String> totalGenes = new HashSet<String>();
        Set<String> totalMiRNAs = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            //System.out.println(line);
            if (!tokens[2].equals("human"))
                continue;
            fis.add(tokens[3] + " " + tokens[7]);
            totalGenes.add(tokens[3]);
            totalMiRNAs.add(tokens[7]);
        }
        fu.close();
        System.out.println("Total interctions: " + fis.size());
        System.out.println("Total genes: " + totalGenes.size());
        System.out.println("Total miRNAs: " + totalMiRNAs.size());
    }
    
    @Test
    public void checkTFs() throws IOException {
        String fileName = "datasets/TransfactionFactors/nrg2538-s3_table.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        int c = 0;
        Set<String> tfs = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens[0].equals("a") ||
                tokens[0].equals("b") ||
                tokens[0].equals("other")) {
                c ++;
                if (tokens.length > 5)
                    tfs.add(tokens[5]);
            }
        }
        fu.close();
        System.out.println("Total lines: " + c);
        System.out.println("Total TFs: " + tfs.size());
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total genes in FI: " + fiGenes.size());
        Set<String> shared = InteractionUtilities.getShared(fiGenes, tfs);
        System.out.println("Shared: " + shared.size());
        fileName = "datasets/TransfactionFactors/TFsInEncode.txt";
        fu.setInput(fileName);
        Set<String> encodeTFs = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            line = line.trim();
            if (line.toLowerCase().contains("input"))
                continue; // Control
            if (line.contains("(")) {
                int index = line.indexOf("(");
                line = line.substring(0, index).trim();
            }
            encodeTFs.add(line);
        }
        System.out.println("Total encode TFs: " + encodeTFs.size());
        shared = InteractionUtilities.getShared(fiGenes, encodeTFs);
        System.out.println("Shared in TFs: " + shared.size());
    }
    
}
