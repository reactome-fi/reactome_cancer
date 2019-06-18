/*
 * Created on Sep 22, 2009
 *
 */
package org.reactome.cancer;

import static org.reactome.cancer.NatureGBMAnalyzer.TCGA_GBM_DIR;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.r3.graph.BreadthFirstSearch;
import org.reactome.r3.graph.BreadthFirstSearch.Edge;
import org.reactome.r3.graph.BreadthFirstSearch.TreeNode;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.MathUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to do TCGA GBM human curated pathway related analysis.
 * @author wgm
 *
 */
public class NatureGBMPathwayAnalyzer {
    private FileUtility fu = new FileUtility();
    private NatureGBMAnalyzer gbmAnalyzer;
    private String sifFileName = TCGA_GBM_DIR + "KnownGBMPathways.sif";
    private String nodeTypeFileName = TCGA_GBM_DIR + "GeneListInGBMPathway.txt";
    
    public NatureGBMPathwayAnalyzer() {
        gbmAnalyzer = new NatureGBMAnalyzer();
    }
    

    @Test
    public void searchForFunctionalExamples() throws Exception {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        Map<String, Set<String>> geneToPartners = new BreadthFirstSearch().generateIdToPartnersMap(fis);
        // Target genes
        Set<String> interactionsInPathway = fu.loadInteractionsInSifFile(sifFileName);
        Set<String> pathwayGenes = InteractionUtilities.grepIDsFromInteractions(interactionsInPathway);
        System.out.println("Total pathway genes: " + pathwayGenes.size());
        // Load altered genes
        List<String> samples = gbmAnalyzer.loadResequencedSamples();
        Set<String> alteredGenes = gbmAnalyzer.getAlteredGenesInSamples(samples);
        System.out.println("Total altered genes: " + alteredGenes.size());
        Set<String> allGenes = new HashSet<String>(pathwayGenes);
        allGenes.addAll(alteredGenes);
        Set<String> allFIs = InteractionUtilities.getFIs(allGenes, fis);
        Set<String> addedGenes = InteractionUtilities.grepIDsFromInteractions(allFIs);
        addedGenes.removeAll(pathwayGenes);
        System.out.println("Total added genes: " + addedGenes.size());
        
        Map<String, Set<String>> sampleToMutations = gbmAnalyzer.loadSampleToMutations();
        Map<String, Set<String>> mutationToSamples = InteractionUtilities.switchKeyValues(sampleToMutations);
        Map<String, Set<String>> sampleToCNVAmpGenes = gbmAnalyzer.loadSampleToCNVAmplifiedGenes();
        Map<String, Set<String>> cnvAmpGenesToSamples = InteractionUtilities.switchKeyValues(sampleToCNVAmpGenes);
        System.out.println("Total amplifed genes: " + cnvAmpGenesToSamples.size());
        Map<String, Set<String>> sampleToCNVDelGenes = gbmAnalyzer.loadSampleToCNVDeletedGenes();
        Map<String, Set<String>> cnvDelGenesToSamples = InteractionUtilities.switchKeyValues(sampleToCNVDelGenes);
        System.out.println("Total deleted genes: " + cnvDelGenesToSamples.size());
        System.out.println("Gene\tMutated_Samples\tAmplified_Samples\tDeleted_Samples\tFI_Pathway_Genes\tList");
        //singleGenes.addAll(checkingGenes);
        for (String gene : addedGenes) {
            Set<String> mutatedSamples = mutationToSamples.get(gene);
            Set<String> cnvAmpSamples = cnvAmpGenesToSamples.get(gene);
            Set<String> cnvDelSamples = cnvDelGenesToSamples.get(gene);
            // Get the total interactions
            Set<String> partners = geneToPartners.get(gene);
            partners.retainAll(pathwayGenes);
            System.out.println(gene + "\t" + 
                               (mutatedSamples == null ? 0 : mutatedSamples.size()) + "\t" +
                               (cnvAmpSamples == null ? 0 : cnvAmpSamples.size()) + "\t" + 
                               (cnvDelSamples == null ? 0 : cnvDelSamples.size()) + "\t" +
                               partners.size() + "\t" + 
                               partners);
                                       
        }
    }
    
    /**
     * This method is used to check how many FIs can be added to the curated GBM pathway.
     * @throws IOException
     */
    @Test
    public void checkAddableFIsToPathway() throws Exception {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        // Target genes
        Set<String> interactionsInPathway = fu.loadInteractionsInSifFile(sifFileName);
        Set<String> pathwayGenes = InteractionUtilities.grepIDsFromInteractions(interactionsInPathway);
        System.out.println("Total pathway genes: " + pathwayGenes.size());
        // Load altered genes
        List<String> samples = gbmAnalyzer.loadResequencedSamples();
        Set<String> alteredGenes = gbmAnalyzer.getAlteredGenesInSamples(samples);
        System.out.println("Total altered genes: " + alteredGenes.size());
        Set<String> pathwayAndAlteredGenes = new HashSet<String>();
        pathwayAndAlteredGenes.addAll(alteredGenes);
        pathwayAndAlteredGenes.addAll(pathwayGenes);
        // Get FIs for these genes
        Set<String> fisInPathwayAndAlteredGenes = InteractionUtilities.getFIs(pathwayAndAlteredGenes, 
                                                                              fis);
        System.out.println("Total FIs for pathway and altered genes: " + fisInPathwayAndAlteredGenes.size());
        // Remove known FIs to get new FIs
        // Get the original FIs
        Set<String> fisInPathway = extractFIsFromGBMPathway();
        System.out.println("Known FIs from pathway: " + fisInPathway.size());
        fisInPathwayAndAlteredGenes.removeAll(fisInPathway);
        System.out.println("New FIs for pathway and altered genes: " + fisInPathwayAndAlteredGenes.size());
        // Want to count FIs between new genes and pathway genes only
        Set<String> fisBetweenNewAndOld = new HashSet<String>();
        alteredGenes.removeAll(pathwayGenes);
        for (String fi : fisInPathwayAndAlteredGenes) {
            String[] partners = fi.split("\t");
            if (alteredGenes.contains(partners[0]) && pathwayGenes.contains(partners[1]))
                fisBetweenNewAndOld.add(fi);
            else if (alteredGenes.contains(partners[1]) && pathwayGenes.contains(partners[0]))
                fisBetweenNewAndOld.add(fi);
        }
        System.out.println("New FIs between pathway genes and newly added FIs: " + fisBetweenNewAndOld.size());
//        // Want to output
//        fu.saveInteractions(fisBetweenNewAndOld,
//                            TCGA_GBM_DIR + "FIsForKnownPathways.txt");
//        // Want to check how many of these genes are predicted
//        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fisBetweenNewAndOld);
//        genes.removeAll(pathwayGenes);
//        System.out.println("Total new genes added: " + genes.size());
//        HibernateFIReader fiReader = new HibernateFIReader();
//        SessionFactory sf = fiReader.initSession();
//        Session session = sf.openSession();
//        Set<String> predicted = new HashSet<String>();
//        for (String fi : fisBetweenNewAndOld) {
//            String[] tokens = fi.split("\t");
//            List<Interaction> fiList = fiReader.queryFIsBasedOnGeneNames(tokens[0], 
//                                                                         tokens[1],
//                                                                         session);
//            if (fiList.size() != 1)
//                continue;
//            Interaction fi1 = fiList.get(0);
//            if (fi1.getEvidence() != null) {
//                predicted.add(fi);
//            }
//        }
//        System.out.println("Total predicted FIs between pathway genes and newly added FIs: " + predicted.size());
//        session.close();
    }
    
    private Map<String, Set<String>> loadNodeTypeToEntities() throws IOException {
        // Counting
        Map<String, Set<String>> typeToEntities = new HashMap<String, Set<String>>();
        fu.setInput(nodeTypeFileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (line.trim().length() == 0)
                break;
            String[] tokens = line.split("\t");
            String gene = tokens[0];
            if (gene.contains(":")) {
                int index = gene.indexOf(":");
                gene = gene.substring(0, index);
            }
            String type = tokens[1];
            Set<String> entities = typeToEntities.get(type);
            if (entities == null) {
                entities = new HashSet<String>();
                typeToEntities.put(type, entities);
            }
            entities.add(gene);
        }
        fu.close();
        return typeToEntities;
    }
    
    /**
     * This method is used to check FIs in the curated pathway.
     * @throws IOException
     */
    @Test
    public void checkFIsInPathway() throws IOException {
        // Get the original PPIs among proteins
        Set<String> fisInPathway = extractFIsFromGBMPathway();
        System.out.println("Curated functional interactions: " + fisInPathway.size());
        // Check how many these FIs are in our FI network
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_FILE_NAME);
        Set<String> copy = new HashSet<String>(fisInPathway);
        copy.retainAll(fis);
        System.out.println("FIs in the FI network: " + copy.size() + " " + (double)copy.size() / fisInPathway.size());
        copy = new HashSet<String>(fisInPathway);
        copy.removeAll(fis);
        System.out.println("FIs not in the FI network: " + copy.size());
        for (String fi : copy)
            System.out.println(fi);
        // Target genes
        Set<String> pathwayGenes = InteractionUtilities.grepIDsFromInteractions(fisInPathway);
        System.out.println("Total pathway genes: " + pathwayGenes.size());
        Set<String> fisInPathwayGenes = InteractionUtilities.getFIs(pathwayGenes, fis);
        System.out.println("Total FIs in pathway genes from the FI network: " + fisInPathwayGenes.size());
        fisInPathwayGenes.removeAll(fisInPathway);
        System.out.println("Total FIs not in the human curated pathway: " + fisInPathwayGenes.size());
    }
    
    /**
     * This method is used to extract FIs from a human curated GBM pathway as defined by our FI.
     * @param relationshipsInPathway
     * @param typeToEntities
     * @return
     * @throws IOException
     */
    private Set<String> extractFIsFromGBMPathway() throws IOException {
        Map<String, Set<String>> typeToEntities = loadNodeTypeToEntities();
        // Load proteins only from the human curated pathways
//        for (String type : typeToEntities.keySet()) {
//            Set<String> entities = typeToEntities.get(type);
//            System.out.println(type + ": " + entities.size());
//        }
        // Need to figure out protein families to members relationships
        Map<String, Set<String>> familyToMembers = new HashMap<String, Set<String>>();
        // Also need to figure protein complex to subunit relationships
        Map<String, Set<String>> complexToSubunits = new HashMap<String, Set<String>>();
        // Original interactions to be parsed: types of activates, inhibits, 
        Set<String> interactions = new HashSet<String>();
        Set<String> interactionTypes = new HashSet<String>();
        // The above two maps should be figured out from the original sif file
        String sifFile = TCGA_GBM_DIR + "KnownGBMPathways.sif";
        fu.setInput(sifFile);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            interactionTypes.add(tokens[1]);
            // Check type
            if (tokens[1].equals("includes")) {
                // Protein family-member relationship
                InteractionUtilities.addElementToSet(familyToMembers, tokens[0], tokens[2]);
            }
            else if (tokens[1].equals("contains")) {
                InteractionUtilities.addElementToSet(complexToSubunits, tokens[0], tokens[2]);
            }
            else if (tokens[1].equals("inhibits") || tokens[1].equals("activates")) {
                interactions.add(tokens[0] + "\t" + tokens[2]);
            }
        }
        fu.close();
        //System.out.println("Interaction types: " + interactionTypes);
        //System.out.println("Total interactions before expanding: " + interactions.size());
        Set<String> fis = new HashSet<String>();
        Set<String> proteins = typeToEntities.get("protein");
        Set<String> complexes = typeToEntities.get("complex");
        Set<String> families = typeToEntities.get("family");
        for (String fi : interactions) {
            int index = fi.indexOf("\t");
            String partner1 = fi.substring(0, index);
            String partner2 = fi.substring(index + 1);
            index = partner1.indexOf(":");
            if (index > 0)
                partner1 = partner1.substring(0, index);
            index = partner2.indexOf(":");
            if (index > 0)
                partner2 = partner2.substring(0, index);
            // Check the types
            List<String> list1 = getListForFIPartner(partner1,
                                                     typeToEntities,
                                                     familyToMembers,
                                                     complexToSubunits);
            List<String> list2 = getListForFIPartner(partner2,
                                                     typeToEntities,
                                                     familyToMembers,
                                                     complexToSubunits);
            Set<String> extractedFIs = generateFIs(list1, list2);
            fis.addAll(extractedFIs);
        }
        return fis;
    }
    
    private Set<String> generateFIs(List<String> list1,
                                    List<String> list2) {
        if (list1.size() == 0 || list2.size() == 0)
            return new HashSet<String>();
        Set<String> fis = new HashSet<String>();
        for (String protein1 : list1) {
            for (String protein2 : list2) {
                fis.add(InteractionUtilities.generateFIFromGene(protein1, protein2));
            }
        }
        return fis;
    }
    
    private List<String> getListForFIPartner(String partner,
                                             Map<String, Set<String>> typeToEntities,
                                             Map<String, Set<String>> familyToMembers,
                                             Map<String, Set<String>> complexToSubunits) {
        Set<String> set = new HashSet<String>();
        getListForPartnerRecursively(partner,
                                     typeToEntities,
                                     familyToMembers,
                                     complexToSubunits,
                                     set);
        return new ArrayList<String>(set);
    }
    
    private void getListForPartnerRecursively(String partner, 
                                              Map<String, Set<String>> typeToEntities,
                                              Map<String, Set<String>> familyToMembers,
                                              Map<String, Set<String>> complexToSubunits,
                                              Set<String> set) {
        // Check the stop condition
        Set<String> proteins = typeToEntities.get("protein");
        if (proteins.contains(partner)) {
            set.add(partner);
            return;
        }
        Set<String> complexes = typeToEntities.get("complex");
        if (complexes.contains(partner)) {
            Set<String> subunits = complexToSubunits.get(partner);
            if (subunits == null) {
                System.err.println(partner + " as no complex subunits!");
            }
            else {
                for (String subunit : subunits)
                    getListForPartnerRecursively(subunit,
                                                 typeToEntities, 
                                                 familyToMembers, 
                                                 complexToSubunits,
                                                 set);
            }
            return;
        }
        Set<String> families = typeToEntities.get("family");
        if (families.contains(partner)) {
            Set<String> members = familyToMembers.get(partner);
            if (members == null)
                System.err.println(partner + " has no family member!");
            else {
                for (String member : members)
                    getListForPartnerRecursively(member,
                                                 typeToEntities,
                                                 familyToMembers, 
                                                 complexToSubunits, 
                                                 set);
            }
        }
    }

    @Test
    public void checkMutationInKnownPathway() throws Exception {
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total FI genes (in biggest component): " + fiGenes.size());
        // Target genes
        Set<String> interactionsInPathway = fu.loadInteractionsInSifFile(sifFileName);
        Set<String> pathwayGenes = InteractionUtilities.grepIDsFromInteractions(interactionsInPathway);
        System.out.println("Total pathway genes: " + pathwayGenes.size());
        // Get altered genes
        List<String> samples = gbmAnalyzer.loadResequencedSamples();
        Set<String> alteredGenes = gbmAnalyzer.getAlteredGenesInSamples(samples);
        // Altered genes in known pathway
        Set<String> copy = new HashSet<String>(pathwayGenes);
        copy.retainAll(alteredGenes);
        double percent = (double) copy.size() / pathwayGenes.size();
        System.out.println("Pathway genes have been altered: " + copy.size() + " (" + percent + ")");
        // Pathway genes in our FI network
        copy = new HashSet<String>(pathwayGenes);
        copy.retainAll(fiGenes);
        System.out.println("Pathway genes in the FI network: " + copy.size());
        int size = copy.size();
        copy.retainAll(alteredGenes);
        percent = (double) copy.size() / size;
        System.out.println("Pathway genes in the FI network and altered: " + copy.size() + " (" + percent + ")");
    }

    @Test
    public void analyzeFIsWithKnownGBMPathway() throws Exception {
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total FI genes (in biggest component): " + fiGenes.size());
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        // Target genes
        Set<String> interactionsInPathway = fu.loadInteractionsInSifFile(sifFileName);
        Set<String> pathwayGenes = InteractionUtilities.grepIDsFromInteractions(interactionsInPathway);
        System.out.println("Total pathway genes: " + pathwayGenes.size());
        Set<String> pathwayGenesCopy = new HashSet<String>(pathwayGenes);
        pathwayGenes.retainAll(fiGenes);
        System.out.println("Pathway genes in FI network: " + pathwayGenes.size());
        // Load altered genes
        List<String> samples = gbmAnalyzer.loadResequencedSamples();
        Set<String> alteredGenes = gbmAnalyzer.getAlteredGenesInSamples(samples);
        System.out.println("Total altered genes: " + alteredGenes.size());
        Set<String> copy = new HashSet<String>(pathwayGenesCopy);
        copy.retainAll(alteredGenes);
        System.out.println("Altered genes in original pathawy: " + copy.size());
        copy.retainAll(fiGenes);
        System.out.println("Altered genes in both original pathway and FI network: " + copy.size());
        alteredGenes.retainAll(fiGenes);
        System.out.println("Altered genes in FI network: " + alteredGenes.size());
        // Don't count genes already known
        alteredGenes.removeAll(pathwayGenes);
        System.out.println("Altered genes in FI network but not in pathway: " + alteredGenes.size());
        // Get known pathway partners
        Set<String> fiPartnersToKnown = new HashSet<String>();
        for (String gene : pathwayGenes)
            fiPartnersToKnown.addAll(geneToPartners.get(gene));
        fiPartnersToKnown.removeAll(pathwayGenes);
        copy = new HashSet<String>(fiPartnersToKnown);
        copy.retainAll(alteredGenes);
        double pvalue = MathUtilities.calculateHypergeometricPValue(fiGenes.size() - pathwayGenes.size(), 
                                                                    fiPartnersToKnown.size(),
                                                                    alteredGenes.size(), 
                                                                    copy.size());
        System.out.println("FI partners to known pathway genes: " + fiPartnersToKnown.size());
        System.out.println("Shared: " + copy.size());
        System.out.println("pvalue from hypergeometric test: " + pvalue);
        // Want to generate these FIs interacting with known pathway components
        Set<String> pathwayFIs = new HashSet<String>();
        // Generate map from external partners to Known genes
        final Map<String, Set<String>> newToOld = new HashMap<String, Set<String>>();
        for (String gene : pathwayGenes) {
            Set<String> partners = geneToPartners.get(gene);
            // Remove pathway genes
            copy = new HashSet<String>(partners);
            copy.removeAll(pathwayGenes);
            copy.retainAll(alteredGenes);
            // Generate FIs
            for (String gene1 : copy) {
                String fi = InteractionUtilities.generateFIFromGene(gene, gene1);
                pathwayFIs.add(fi);
                InteractionUtilities.addElementToSet(newToOld,
                                                     gene1,
                                                     gene);
            }
        }
        //String outFileName = TCGA_GBM_DIR + "FIsForKnownPathways.txt";
        //fu.saveInteractions(pathwayFIs, outFileName);
        // Check the size
        List<String> newGenes = new ArrayList<String>(newToOld.keySet());
        System.out.println("Total new genes added to the GBM pathway: " + newGenes.size());
        Collections.sort(newGenes, new Comparator<String>() {
            public int compare(String gene1, String gene2) {
                Set<String> set1 = newToOld.get(gene1);
                Set<String> set2 = newToOld.get(gene2);
                return set2.size() - set1.size();
            }
        });
        for (String gene : newGenes) {
            Set<String> partners = geneToPartners.get(gene);
            Set<String> old = newToOld.get(gene);
            // calculate p-value
            pvalue = MathUtilities.calculateHypergeometricPValue(fiGenes.size(), 
                                                                 partners.size(),
                                                                 pathwayGenes.size(),
                                                                 old.size());
            if (pvalue < 0.0)
                pvalue = 0.0;
            System.out.println(gene + "\t" + old.size() + "\t" + partners.size() + "\t" + pvalue);
        }
        // Want to get the first ten and last ten genes for analysis
        Set<String> newFIs = new HashSet<String>();
        Set<String> checkingGenes = new HashSet<String>();
        // Get the top 10
        for (int i = 0; i < 10; i ++) {
            String newGene = newGenes.get(i);
            Set<String> set = newToOld.get(newGene);
            for (String old : set) {
                String fi = InteractionUtilities.generateFIFromGene(old, newGene);
                newFIs.add(fi);
            }
            checkingGenes.add(newGene);
        }
        // Get the last 10
        for (int i = newGenes.size() - 1; i > newGenes.size() - 11; i--) {
            String newGene = newGenes.get(i);
            Set<String> set = newToOld.get(newGene);
            for (String old : set) {
                String fi = InteractionUtilities.generateFIFromGene(old, newGene);
                newFIs.add(fi);
            }
            checkingGenes.add(newGene);
        }
        //String outFileName = TCGA_GBM_DIR + "FIsForKnownPathwaysTopBottom10.txt";
        //fu.saveInteractions(newFIs, outFileName);
        // In this file, only single FI will be printed out to simplify the view
        newFIs.clear();
        int c = 0;
        Set<String> singleGenes = new HashSet<String>();
        for (String gene : newGenes) {
            Set<String> set = newToOld.get(gene);
            if (set.size() > 1) {
                continue;
            }
            String old = set.iterator().next();
            c ++;
            String fi = InteractionUtilities.generateFIFromGene(old, gene);
            newFIs.add(fi);
            singleGenes.add(gene);
        }
        System.out.println("Total single interacting genes: " + c);
        //String outFileName = TCGA_GBM_DIR + "FIsForKnownPathwaySingle.txt";
        //fu.saveInteractions(newFIs, outFileName);
        // Want to generate FIs from predictions only
//        newFIs.clear();
//        HibernateFIReader fiReader = new HibernateFIReader();
//        SessionFactory sf = fiReader.initSession();
//        Session session = sf.openSession();
//        for (String gene : newGenes) {
//            Set<String> set = newToOld.get(gene);
//            //if (set.size() > 1)
//            //    continue;
//            for (String old : set) {
//                List<Interaction> fiList = fiReader.queryFIsBasedOnGeneNames(gene, 
//                                                                             old,
//                                                                             session);
//                if (fiList.size() != 1)
//                    continue;
//                Interaction fi = fiList.get(0);
//                if (fi.getEvidence() != null) {
//                    newFIs.add(InteractionUtilities.generateFIFromGene(old, gene));
//                }
//            }
//        }
//        System.out.println("Total predicted FIs: " + newFIs.size());
//        session.close();
        //String outFileName = TCGA_GBM_DIR + "FIsForKnownPathwayPredicted.txt";
        //fu.saveInteractions(newFIs, outFileName);
        // Check types for single genes
        Map<String, Set<String>> sampleToMutations = gbmAnalyzer.loadSampleToMutations();
        Map<String, Set<String>> mutationToSamples = InteractionUtilities.switchKeyValues(sampleToMutations);
        Map<String, Set<String>> sampleToCNVAmpGenes = gbmAnalyzer.loadSampleToCNVAmplifiedGenes();
        Map<String, Set<String>> cnvAmpGenesToSamples = InteractionUtilities.switchKeyValues(sampleToCNVAmpGenes);
        Map<String, Set<String>> sampleToCNVDelGenes = gbmAnalyzer.loadSampleToCNVDeletedGenes();
        Map<String, Set<String>> cnvDelGenesToSamples = InteractionUtilities.switchKeyValues(sampleToCNVDelGenes);
        System.out.println("Gene\tMutated_Samples\tAmplified_Samples\tDeleted_Samples");
        //singleGenes.addAll(checkingGenes);
        for (String gene : singleGenes) {
            Set<String> mutatedSamples = mutationToSamples.get(gene);
            Set<String> cnvAmpSamples = cnvAmpGenesToSamples.get(gene);
            Set<String> cnvDelSamples = cnvDelGenesToSamples.get(gene);
            System.out.println(gene + "\t" + 
                               (mutatedSamples == null ? 0 : mutatedSamples.size()) + "\t" +
                               (cnvAmpSamples == null ? 0 : cnvAmpSamples.size()) + "\t" + 
                               (cnvDelSamples == null ? 0 : cnvDelSamples.size()));
                                       
        }
    }
    
    @Test
    public void checkTopologicalPropertiesFromAlteredGenesToKnownGBMPathway() throws Exception {
        // Set up graph
        boolean useShortestPath = true;
        BreadthFirstSearch bfs = new BreadthFirstSearch();
        Set<String> fis = fu.loadInteractions(R3Constants.GENE_FI_BIG_COMP_FILE_NAME);
        Set<String> fiGenes = InteractionUtilities.grepIDsFromInteractions(fis);
        System.out.println("Total FI genes (in biggest component): " + fiGenes.size());
        Map<String, Set<String>> geneToPartners = bfs.generateIdToPartnersMap(fis);
        Map<TreeNode, List<Edge>> nodeToEdges = bfs.initGraph(fis);
        // Target genes
        String sifFileName = TCGA_GBM_DIR + "KnownGBMPathways.sif";
        Set<String> interactionsInPathway = fu.loadInteractionsInSifFile(sifFileName);
        Set<String> pathwayGenes = InteractionUtilities.grepIDsFromInteractions(interactionsInPathway);
        pathwayGenes.retainAll(fiGenes);
        System.out.println("Pathway genes in FI network: " + pathwayGenes.size());
        // Calculate average path among this set of genes
        double avgPath = gbmAnalyzer.calculateShortestPath(new ArrayList<String>(pathwayGenes),
                                                           bfs, 
                                                           nodeToEdges);
        System.out.println("Shortest path among pathway genes: " + avgPath);
        // Load altered genes
        List<String> samples = gbmAnalyzer.loadResequencedSamples();
        Set<String> alteredGenes = gbmAnalyzer.getAlteredGenesInSamples(samples);
        alteredGenes.retainAll(fiGenes);
        System.out.println("Altered genes in FI network: " + alteredGenes.size());
        // Don't count genes already known
        alteredGenes.removeAll(pathwayGenes);
        System.out.println("Altered genes in FI network but not in pathway: " + alteredGenes.size());
        // Start calculation
        int totalPathLength = 0;
        int totalPath = 0;
        if (useShortestPath) {
            // Want to pick up the shortest path between the checking gene and all target genes
            for (String gene : alteredGenes) {
                Map<String, Integer> geneToPath = bfs.getDistances(gene, 
                                                                   pathwayGenes,
                                                                   geneToPartners);
                int currentPath = Integer.MAX_VALUE;
                for (String target : geneToPath.keySet()) {
                    Integer path = geneToPath.get(target);
                    if (path < currentPath)
                        currentPath = path;
                }
                //System.out.println("path for " + gene + ": " + currentPath);
                totalPathLength += currentPath;
                totalPath ++;
            }
        }
        else {
            // Want to pick up all shortest paths between the checking gene and all target genes
            for (String gene : alteredGenes) {
                Map<String, Integer> geneToPath = bfs.getDistances(gene, 
                                                                   pathwayGenes,
                                                                   geneToPartners);
                for (String target : geneToPath.keySet()) {
                    Integer path = geneToPath.get(target);
                    System.out.println(gene + " " + target +": " + path);
                    totalPath ++;
                    totalPathLength += path;
                }
            }
        }
        System.out.println("Total path: " + totalPath);
        System.out.println("Average shortest path: " + totalPathLength / (double)totalPath);
        System.out.println("Total FI genes: " + fiGenes.size());     
        // For permutation
        int permutation = 100;
        List<Double> list = new ArrayList<Double>();
        // For genes
        //      Set<String> sampleGenes = new HashSet<String>(fiGenes);
        //      sampleGenes.removeAll(pathwayGenes);
        // Get 601 sequenced genes
        Set<String> sampleGenes = gbmAnalyzer.loadAllNatureGBM601Genes();
        sampleGenes.retainAll(fiGenes);
        System.out.println("601 genes in FI: " + sampleGenes.size());
        sampleGenes.removeAll(pathwayGenes);
        System.out.println("Genes for sampling: " + sampleGenes.size());
        for (int i = 0; i < permutation; i++) {
            Set<String> permutationGenes = MathUtilities.randomSampling(sampleGenes, alteredGenes.size());
            totalPath = 0;
            totalPathLength = 0;
            // Want to pick up the shortest path between the checking gene and all target genes
            for (String gene : permutationGenes) {
                Map<String, Integer> geneToPath = bfs.getDistances(gene, 
                                                                   pathwayGenes,
                                                                   geneToPartners);
                int currentPath = Integer.MAX_VALUE;
                for (String target : geneToPath.keySet()) {
                    Integer path = geneToPath.get(target);
                    if (path < currentPath)
                        currentPath = path;
                }
                //System.out.println("path for " + gene + ": " + currentPath);
                totalPathLength += currentPath;
                totalPath ++;
            }
            list.add((double)totalPathLength / totalPath);
        }
        Collections.sort(list);
        System.out.println("\nPermutation results:");
        for (int i = 0; i < list.size(); i++)
            System.out.println((i + 1) + "\t" + list.get(i));
    }

}
