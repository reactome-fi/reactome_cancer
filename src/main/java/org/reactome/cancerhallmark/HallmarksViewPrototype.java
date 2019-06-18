/*
 * Created on May 2, 2016
 *
 */
package org.reactome.cancerhallmark;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeModel;

import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.util.SwingImageCreator;
import org.junit.Test;
import org.reactome.r3.graph.GraphAnalyzer;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;

/**
 * @author gwu
 *
 */
public class HallmarksViewPrototype {
    
    public HallmarksViewPrototype() {
    }
    
    @Test
    public void createTabbedPanes() throws Exception {
        JFrame frame = new JFrame("Tabbed Pane");
        JTabbedPane tabbedPane = new JTabbedPane();
        tabbedPane.setTabPlacement(JTabbedPane.BOTTOM);
        frame.getContentPane().add(tabbedPane, BorderLayout.CENTER);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        String[] titles = new String[] {
                "Pathway Overview",
                "Pathway Diagram",
                "Data View",
                "Reaction Network"
        };
        for (String title : titles) {
            tabbedPane.add(title, new JPanel());
        }
        frame.setSize(500, 400);
        frame.setVisible(true);
        
        Thread.sleep(10000000);
    }
    
    @Test
    public void generateHallmarksTree() throws Exception {
        String[] hallmarks = new String[] {
                "Activating invasion and metastasis (AIM)",
                "Evading growth suppressors (EGS)",
                "Evading immune destruction (EIM)",
                "Enabling replicative immortality (ERI)",
                "Genome Instability and Mutation (GIM)",
                "Inducing angiogenesis (IAG)",
                "Resisting cell death (RCD)",
                "Reprogramming of energy metabolism (REM)",
                "Sustaining proliferative Signaling (SPS)",
                "Tumor-promoting Inflammation (TPI)"
        };
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_55_plus_i",
                                            "root",
                                            "macmysql01");
        Collection<GKInstance> c = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Pathway, 
                                                                ReactomeJavaConstants._displayName,
                                                                "=",
                                                                "Programmed Cell Death");
        GKInstance inst = c.iterator().next();
        List<GKInstance> subpathways = inst.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
        for (GKInstance subpathway : subpathways)
            System.out.println(subpathway);
        JTree tree = new JTree();
        tree.setRootVisible(false);
        tree.setShowsRootHandles(true);
        tree.setCellRenderer(new DefaultTreeCellRenderer() {
            private Icon pathwayIcon = new ImageIcon("images/Pathway.gif");
            private Icon reactionIcon = new ImageIcon("images/Reaction.gif");
            
            @Override
            public Component getTreeCellRendererComponent(JTree tree, 
                                                          Object value, 
                                                          boolean selected, 
                                                          boolean expanded,
                                                          boolean leaf, 
                                                          int row, 
                                                          boolean hasFocus) {
                Component comp = super.getTreeCellRendererComponent(tree, value, selected, expanded, leaf, row, hasFocus);
                if (leaf)
                    setIcon(reactionIcon);
                else
                    setIcon(pathwayIcon);
                return comp;
            }
        });
        
        DefaultMutableTreeNode root = new DefaultMutableTreeNode("");
        for (String hallmark : hallmarks) {
            DefaultMutableTreeNode treeNode = new DefaultMutableTreeNode(hallmark);
            root.add(treeNode);
            if (hallmark.equals("Resisting cell death (RCD)")) {
                addTreeNodes(inst, treeNode);
            }
            else
                treeNode.add(new DefaultMutableTreeNode("test")); // Just to show a root
        }
        
        TreeModel treeModel = new DefaultTreeModel(root);
        tree.setModel(treeModel);
        
        JFrame frame = new JFrame();
        frame.getContentPane().add(new JScrollPane(tree), BorderLayout.CENTER);
        frame.setSize(400, 800);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
        Thread.sleep(60* 60 * 1000);
    }
    
    private void addTreeNodes(GKInstance pathway,
                              DefaultMutableTreeNode treeNode) throws Exception {
        if (pathway.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasEvent)) {
            List<GKInstance> hasEvent = pathway.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
            if (hasEvent == null || hasEvent.size() == 0)
                return;
            for (GKInstance subpathway : hasEvent) {
                String name = subpathway.getDisplayName();
                DefaultMutableTreeNode subNode = new DefaultMutableTreeNode(name);
                treeNode.add(subNode);
                addTreeNodes(subpathway, subNode);
            }
        }
    }
    
    @Test
    public void generateCancerHallmarksBars() throws IOException {
        int w = 1000;
        int h = 50;
        w = 50;
        h = 100;
        HallmarksPanel panel = new HallmarksPanel();
        panel.setSize(w, h);
        String fileName = "/Users/Gwu/Dropbox/OHSU/PAR-15-332_CancerHallmarks/Figures/HallmarksBars.pdf";
        fileName = "/Users/Gwu/Dropbox/OHSU/PAR-15-332_CancerHallmarks/Figures/HallmarksBars_Legend.pdf";
        SwingImageCreator.exportImageInPDF(panel,
                                           new File(fileName));
    }
    
    @Test
    public void generateDataPanel() throws IOException {
        int w = 2000;
        int h = 500;
        DatasetPanel panel = new DatasetPanel();
        panel.setSize(w, h);
        String fileName = "/Users/Gwu/Dropbox/OHSU/PAR-15-332_CancerHallmarks/Figures/DataPanel_Focused.pdf";
        SwingImageCreator.exportImageInPDF(panel,
                                           new File(fileName));
    }
    
    private class DatasetPanel extends JComponent {
        private Color posColor = Color.RED;
        private Color negColor = Color.WHITE;
        
        @Override
        public void paint(Graphics g) {
            int w = getWidth();
            int h = getHeight();
            int wcell = 1000;
            int hcell = 100;
            // Set of frequently mutated genes
            Set<Integer> focusedPoints = new HashSet<Integer>();
            focusedPoints.add(250);
            focusedPoints.add(275);
            focusedPoints.add(260);
            Graphics2D g2 = (Graphics2D) g;
            g2.clearRect(0, 0, w, h);
            // Divide into 100 * 500 cells
            double ws = (double) w / wcell;
            double hs = (double) h / hcell;
            double x = 0.0d, y = 0.0d;
            for (int i = 0; i < wcell; i++) {
                for (int j = 0; j < hcell; j++) {
                    Rectangle2D rect = new Rectangle2D.Double(x, y, ws, hs);
                    if (Math.random() > 0.975d ||
                        (focusedPoints.contains(i) && Math.random() > 0.25d)) {
                        g2.setPaint(posColor);
                        g2.fill(rect);
                    }
                    // Don't do anything for negative
                    y += hs;
                }
                y = 0.0d;
                x += ws;
            }
        }
            
    }
    
    private class HallmarksPanel extends JComponent {
        private String colorText = "0,0,153;51,153,0;153,0,153;102,102,0;0,102,102;204,153,0;204,0,102;204,204,255;0,102,51;0,204,204";
        private List<Color> colors;
        
        public HallmarksPanel() {
            colors = new ArrayList<Color>();
            String[] tokens = colorText.split(";");
            for (String token : tokens) {
                String[] tokens1 = token.split(",");
                int r = new Integer(tokens1[0]);
                int g = new Integer(tokens1[1]);
                int b = new Integer(tokens1[2]);
                colors.add(new Color(r, g, b));
            }
        }
        
        @Override
        public void paint(Graphics g) {
            int buffer = 2;
            int w = getWidth();
            int h = getHeight();
            Graphics2D g2 = (Graphics2D) g;
            BasicStroke stroke = new BasicStroke(2.0f);
            g2.setStroke(stroke);
            g2.clearRect(0, 0, w, h);
            int total = 11;
            double step = (double) h / total;
            double y = 0;
            BasicStroke whiteStroke = new BasicStroke(3.0f);
            for (int i = 0; i < total - 1; i++) {
                y += step;
                Line2D line = new Line2D.Double(buffer, y, w - buffer, y);
                g2.setPaint(colors.get(i));
                g2.draw(line);
                if (true)
                    continue;
                // Create randomly broken to simulate 
                int broken = (int) (20 * Math.random()) + 10;
                g2.setPaint(Color.white);
                g2.setStroke(whiteStroke);
                for (int j = 0; j < broken; j++) {
                    int x1 = (int) (w * Math.random());
                    int w1 = (int) (w / 8 * Math.random());
                    line = new Line2D.Double(x1 + buffer, y, x1 + w1 - buffer, y);
                    g2.draw(line);
                }
                g2.setStroke(stroke);
            }
        }

    }
    
    private List<String> selectGenesRandomly(Set<String> fis) {
        GraphAnalyzer graphAnalyzer = new GraphAnalyzer();
        List<Set<String>> components = graphAnalyzer.calculateGraphComponents(fis);
        // Get the largest component
        List<String> geneList = new ArrayList<String>(components.get(0));
        System.out.println("Total genes in the largest component: " + geneList.size());

        List<String> newList = new ArrayList<String>();
        // Extract 1000 genes for quick prototyping
        RandomData randomData = new RandomDataImpl();
        int size = 1750;
//      size = 500;
        int[] indices = randomData.nextPermutation(geneList.size(), size);
        for (int index : indices)
            newList.add(geneList.get(index));
        
        fis = InteractionUtilities.getFIs(newList, fis);
        components = graphAnalyzer.calculateGraphComponents(fis);
        Set<String> rtn = components.get(0);
        System.out.println("Selected genes: " + rtn.size());
        return new ArrayList<String>(rtn);
    }
    
    @Test
    public void generateFIMatrix() throws Exception {
        String dir = "/Users/Gwu/Documents/temp/";
        String src = dir + "FIsInGene_031516_with_annotations.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(src);
        Set<String> fis = new HashSet<String>();
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            line = fu.readLine();
            String[] tokens = line.split("\t");
            fis.add(tokens[0] + "\t" + tokens[1]);
        }
        fu.close();
        System.out.println(fis.size());
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        List<String> geneList = new ArrayList<String>(genes);
        Collections.sort(geneList);
        System.out.println("Total genes: " + geneList.size());
        
        String target = dir + "FIMatrix.txt";
        
        geneList = selectGenesRandomly(fis);
        target = dir + "FIMatrix_1000.txt";
        
        fu.setOutput(target);
        StringBuilder builder = new StringBuilder();
        builder.append("Gene");
        for (String gene : geneList)
            builder.append("\t").append(gene);
        fu.printLine(builder.toString());
        Map<String, Set<String>> geneToPartners = InteractionUtilities.generateProteinToPartners(fis);
        for (String gene : geneList) {
            builder.setLength(0);
            builder.append(gene);
            Set<String> partners = geneToPartners.get(gene);
            for (String gene1 : geneList) {
                builder.append("\t");
                if (partners.contains(gene1))
                    builder.append("1");
                else
                    builder.append("0");
            }
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
}
