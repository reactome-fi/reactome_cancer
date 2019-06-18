/*
 * Created on Jan 27, 2010
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.junit.Test;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * TODO: This class should be copied to Windows machine for OHSU soon.
 * @author wgm
 *
 */
public class ADAnalyzer extends CancerResequenceDataSetAnalyzer {
    private final String DIR_NAME = "/Users/wgm/Documents/wgm/OHSU/";
    private FileUtility fu = new FileUtility();
    
    public ADAnalyzer() {
    }
    
    /**
     * This method is used to annotate genes based on p-values.
     * @throws Exception
     */
    @Test
    public void annotateGenesFromArrayViaPValueCutoff() throws Exception {
        double pvalueCutoff = 0.001d; 
        //double pvalueCutoff = 0.00334433; // FDR = 0.05
        //double pvalueCutoff = 0.030675; // FDR = 0.10
        //double foldChangeCutoff = 1.5;
        double foldChangeCutoff = 0.0d;
        String dirName = DIR_NAME + "pmid19937809/";
        String fileName = dirName + "sm001.txt";
        fu.setInput(fileName);
        Set<String> upGenes = new HashSet<String>();
        Set<String> downGenes = new HashSet<String>();
        String line = fu.readLine();
        Set<String> geneSet = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Double pvalue = new Double(tokens[4]);
            if (pvalue > pvalueCutoff)
                continue;
            geneSet.clear();
            String[] geneTokens = tokens[2].split("( \\/\\/\\/ )");
            // Check fold change
            Double foldChange = new Double(tokens[3]);
            if (foldChange >= foldChangeCutoff) {
                // Up genes
                for (String gene : geneTokens)
                    upGenes.add(gene);
            }
            else if (foldChange <= foldChangeCutoff) {
                // Down genes
                for (String gene : geneTokens)
                    downGenes.add(gene);
            }
        }
        fu.close();
        // Just a check
        System.out.println("P-value cutoff: " + pvalueCutoff);
        System.out.println("Up genes: " + upGenes.size());
        PathwayBasedAnnotator annotator = initAnnotator();
        annotator.annotateGenesWithFDR(upGenes, AnnotationType.Pathway);
        System.out.println("\nDown genes: " + downGenes.size());
        annotator.annotateGenesWithFDR(downGenes, AnnotationType.Pathway);
    }
    
    /**
     * This method is used to annotate genes with pathways.
     * @throws IOException
     */
    @Test
    public void annotateGenesViaPValueCutoff() throws Exception {
        String dirName = DIR_NAME  + "results/";
        String fileNames[] = new String[] {
                "GeneToPValueChisq_0.txt",
                "GeneToPValueChisq_500.txt",
                "GeneToPValueChisq_10000.txt",
                "GeneToPValueChisq_500000.txt",
                "GeneToPValueLogistic_0.txt",
                "GeneToPValueLogistic_500.txt",
                "GeneToPValueLogistic_10000.txt",
                "GeneToPValueLogistic_500000.txt"
        };
        double pvalueCutoff = 0.01;
        PathwayBasedAnnotator annotator = initAnnotator();
        FileUtility fu = new FileUtility();
        System.out.println("P-value cutoff: " + pvalueCutoff);
        for (String fileName : fileNames) {
            String tmp = dirName + fileName;
            fu.setInput(tmp);
            String line = fu.readLine();
            Set<String> genes = new HashSet<String>();
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                Double pvalue = new Double(tokens[1]);
                if (pvalue <= pvalueCutoff)
                    genes.add(tokens[0]);
            }
            fu.close();
            System.out.println("\n" + fileName + ": " + genes.size());
            annotator.annotateGenesWithFDR(genes, AnnotationType.Pathway);
        }
    }
    
    private PathwayBasedAnnotator initAnnotator() {
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setFDRThreshold(0.20);
        annotator.setPValueThreshold(0.05);
        return annotator;
    }
    
    private Set<String> loadADGenes() throws IOException {
        String fileName = DIR_NAME + "topGenesfrom AlzGene.txt";
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> genes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            genes.add(tokens[1]);
        }
        fu.close();
        return genes;
    }
    
    @Test
    public void checkFIsInGenes() throws Exception {
        Set<String> genes = loadADGenes();
        // Generate a list of genes
        StringBuilder builder = new StringBuilder();
        for (String gene : genes)
            builder.append(gene).append(", ");
        System.out.println(builder.toString());
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> shared = InteractionUtilities.getShared(genes, fiGenes);
        Set<String> fisInGenes = InteractionUtilities.getFIs(genes, fis);
        System.out.println("Total genes: " + genes.size());
        System.out.println("Total genes in the FI network: " + shared.size());
        System.out.println("Total FIs for genes: " + fisInGenes.size());
        // Want to create a sub-network for AD genes in the FI network
        String outFIFileName = DIR_NAME + "FIsForTopADGenes.txt";
        String outSpanFileName = DIR_NAME + "SpanForTopADGenes.txt";
        calculateMinimumSpan(outFIFileName, 
                             outSpanFileName,
                             new ArrayList<String>(shared));
    }
    
    @Test
    public void annotateGenes() throws Exception {
        Set<String> genes = loadADGenes();
        System.out.println("Total genes: " + genes.size());
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setFDRThreshold(0.20);
        annotator.setPValueThreshold(0.05);
        System.out.println("\nAnnotate with pathways:");
        annotator.annotateGenesWithFDR(genes, AnnotationType.Pathway);
        System.out.println("\nAnnotate with GO BP:");
        annotator.annotateGenesWithFDR(genes, AnnotationType.BP);
        System.out.println("\nAnnotate with GO MF:");
        annotator.annotateGenesWithFDR(genes, AnnotationType.MF);
        System.out.println("\nAnnotate with GO CC:");
        annotator.annotateGenesWithFDR(genes, AnnotationType.CC);
    }
    
}
