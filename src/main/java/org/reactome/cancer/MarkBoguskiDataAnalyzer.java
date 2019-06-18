/*
 * Created on Nov 5, 2012
 *
 */
package org.reactome.cancer;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Test;
import org.reactome.annotate.AnnotationType;
import org.reactome.annotate.GeneSetAnnotation;
import org.reactome.annotate.PathwayBasedAnnotator;
import org.reactome.r3.util.FileUtility;

/**
 * This class is used to process data from Mark Bogustki.
 * @author gwu
 *
 */
public class MarkBoguskiDataAnalyzer {
    private final String DIR_NAME = "/Users/gwu/Documents/Lincoln/MarkBoguski/";
    
    public MarkBoguskiDataAnalyzer() {
        
    }
    
    @Test
    public void checkExomeMutationFile() throws Exception {
        String fileName = DIR_NAME + "Exome sequence variants.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        Set<String> samples = new HashSet<String>();
        int consequenceIndex = 0;
        int geneIndex = 0;
        int mutPctIndex = 0;
        String[] headers = line.split("\t");
        for (int i = 0; i < headers.length; i++) {
            if (headers[i].equals("Consequence")) {
                consequenceIndex = i;
            }
            else if (headers[i].equals("GeneName")) {
                geneIndex = i;
            }
            else if (headers[i].equals("MutPct")) {
                mutPctIndex = i;
            }
        }
        List<String> pickedConsequences = getPickedConsequences();
        Set<String> consequences = new HashSet<String>();
        Set<String> genes = new HashSet<String>();
        Set<String> homoGenes = new HashSet<String>(); // Actually for mutpct >= 50%
        while ((line = fu.readLine()) != null) {
//            System.out.println(line);
            String[] tokens = line.split("\t");
            samples.add(tokens[0]);
            String consequence = tokens[consequenceIndex];
            int index = consequence.indexOf(":");
            if (index > 0)
                consequence = consequence.substring(0, index);
            consequences.add(consequence);
            if (pickedConsequences.contains(consequence)) {
                genes.add(tokens[geneIndex]);
                // Get rid of last %
                String tmp = tokens[mutPctIndex];
                double mutPct = new Double(tmp.substring(0, tmp.length() - 1)) / 100.0d;
                if (mutPct >= 0.90d)
                    homoGenes.add(tokens[geneIndex]);
            }
        }
        fu.close();
        System.out.println("Total samples: " + samples);
        System.out.println("Total consequences: " + consequences.size());
        for (String consequence : consequences) {
            System.out.println(consequence);
        }
        System.out.println("Total genes: " + genes.size());
        System.out.println("Total homo genes: " + homoGenes.size());
//        for (String gene : homoGenes)
//            System.out.println(gene);
        fileName = DIR_NAME + "NonSynmousGenes.txt";
        fu.saveInteractions(genes, fileName);
        fileName = DIR_NAME + "HomoGenes.txt";
        fu.saveInteractions(homoGenes, fileName);
        PathwayBasedAnnotator annotator = new PathwayBasedAnnotator();
        annotator.setFDRThreshold(0.05d);
//        annotator.setPValueThreshold(0.05d);
        System.out.println("\nAnnotate all genes:");
        List<GeneSetAnnotation> annotations = annotator.annotateGenesWithFDR(genes, AnnotationType.Pathway);
        outputAnnotations(annotations);
        System.out.println("\nAnnotate homo genes (mut > 90%):");
        annotations = annotator.annotateGenesWithFDR(homoGenes, AnnotationType.Pathway);
        outputAnnotations(annotations);
    }
    
    private void outputAnnotations(List<GeneSetAnnotation> annotations) {
        System.out.println("Pathway\tHitNumber\tFDR\tGenes");
        for (GeneSetAnnotation annotation : annotations) {
            if (annotation.getFdr().equals("1.00e+00"))
                continue;
            if (annotation.getFdr().startsWith("<"))
                System.out.println(annotation.getTopic() + "\t" + 
                                   annotation.getHitNumber() + "\t" + 
                                   annotation.getFdr() + "\t" + 
                                   annotation.getHitIds());
            else {
                Double value = new Double(annotation.getFdr());
                if (value <= 0.05)
                    System.out.println(annotation.getTopic() + "\t" +
                                       annotation.getHitNumber() + "\t" + 
                                       annotation.getFdr() + "\t" + 
                                       annotation.getHitIds());
            }
        }
    }
    
    private List<String> getPickedConsequences() {
        String[] consequences = new String[] {
//                "UPSTREAM",
//                "DOWNSTREAM_UTR",
                "START_LOST",
//                "SYNONYMOUS_CODING",
                "NON_SYNONYMOUS_CODING",
//                "DOWNSTREAM",
                "CODON_DELETION",
                "CODON_CHANGE_PLUS_CODON_INSERTION",
//                "UPSTREAM_UTR",
                "STOP_LOST",
                "NON_SYNONYMOUS_START",
                "CODON_INSERTION",
//                "UTR_3_PRIME",
                "SPLICE_SITE_DONOR",
                "STOP_GAINED",
                "SPLICE_SITE_ACCEPTOR",
                "CODON_CHANGE_PLUS_CODON_DELETION",
//                "SYNONYMOUS_STOP",
                "FRAME_SHIFT"
        };
        return Arrays.asList(consequences);
    }
    
}
