/*
 * Created on Mar 7, 2007
 *
 */
package org.reactome.r3;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;
import org.junit.Test;

/**
 * This utility class is used to merge reaction instances.
 * @author guanming
 *
 */
public class ReactionInstanceMerger {
    
    public ReactionInstanceMerger() {
    }
    
    @SuppressWarnings("unchecked")
    public List<GKInstance> merge(List<GKInstance> reactions) throws Exception {
        // Get all participants
        Set<GKInstance> entities = new HashSet<GKInstance>();
        for (GKInstance rxt : reactions) {
            entities.addAll(InstanceUtilities.getReactionParticipants(rxt));
        }
        // Generate keys for entities
        Map<GKInstance, String> entityToKey = new HashMap<GKInstance, String>();
        for (GKInstance entity : entities)
            generateEntityKey(entity, entityToKey);
        Map<String, List<GKInstance>> keyToReactions = new HashMap<String, List<GKInstance>>();
        for (GKInstance rxt : reactions) {
            String key = generateReactionKey(rxt, entityToKey);
            List<GKInstance> list = keyToReactions.get(key);
            if (list == null) {
                list = new ArrayList<GKInstance>();
                keyToReactions.put(key, list);
            }
            list.add(rxt);
        }
        List<GKInstance> rtn = new ArrayList<GKInstance>();
        // Need to pick
        //List<String> keys = new ArrayList<String>(keyToReactions.keySet());
        //Collections.sort(keys);
        //FileUtility fu = new FileUtility();
        //fu.setOutput("results/v2/MergingResultsForTwoApoptosis.txt");
        //for (String key : keys) {
        for (Iterator<String> keys = keyToReactions.keySet().iterator(); keys.hasNext();) {
            String key = keys.next();
            //fu.printLine(key);
            List<GKInstance> list = keyToReactions.get(key);
            if (list.size() == 1)
                rtn.add(list.get(0));
            else {
                // Give reactions from Reactome have higher prority
                GKInstance found = null;
                for (GKInstance rxt : list) {
                    if (rxt.getAttributeValue(ReactomeJavaConstants.dataSource) == null) {
                        found = rxt;
                        break;
                    }
                }
                if (found != null)
                    rtn.add(found);
                else
                    rtn.add(list.get(0));
            }
        }
        //fu.close();
        return rtn;
    }
    
    private String generateReactionKey(GKInstance reaction,
                                       Map<GKInstance, String> entityToKey) throws Exception {
        StringBuilder builder = new StringBuilder();
        List inputs = reaction.getAttributeValuesList(ReactomeJavaConstants.input);
        if (inputs != null && inputs.size() > 0) {
            builder.append("input:");
            for (Iterator it = inputs.iterator(); it.hasNext();) {
                GKInstance input = (GKInstance) it.next();
                String key = entityToKey.get(input);
                builder.append(key);
                if (it.hasNext())
                    builder.append("+");
            }
        }
        List outputs = reaction.getAttributeValuesList(ReactomeJavaConstants.output);
        if (outputs != null && outputs.size() > 0) {
            builder.append("||output:");
            for (Iterator it = outputs.iterator(); it.hasNext();) {
                GKInstance output = (GKInstance) it.next();
                String key = entityToKey.get(output);
                builder.append(key);
                if (it.hasNext())
                    builder.append("+");
            }
        }
        List cas = reaction.getAttributeValuesList(ReactomeJavaConstants.catalystActivity);
        if (cas != null && cas.size() > 0) {
            builder.append("||catalyst:");
            for (Iterator it = cas.iterator(); it.hasNext();) {
                GKInstance ca = (GKInstance) it.next();
                GKInstance catalyst = (GKInstance) ca.getAttributeValue(ReactomeJavaConstants.physicalEntity);
                String key = entityToKey.get(catalyst);
                builder.append(key);
                if (it.hasNext())
                    builder.append("+");
            }
        }
        Collection regulations = reaction.getReferers(ReactomeJavaConstants.regulatedEntity);
        if (regulations != null && regulations.size() > 0) {
            builder.append("||regulator:");
            for (Iterator it = regulations.iterator(); it.hasNext();) {
                GKInstance regulation = (GKInstance) it.next();
                GKInstance regulator = (GKInstance) regulation.getAttributeValue(ReactomeJavaConstants.regulator);
                if (regulator == null) {
                    builder.append("unknown");
                    if (it.hasNext())
                        builder.append("+");
                }
                else if (regulator.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity)) {
                    String key = entityToKey.get(regulator);
                    builder.append(key);
                    if (it.hasNext())
                        builder.append("+");
                }
            }
        }
        return builder.toString();
    }
    
    private void generateEntityKey(GKInstance entity,
                                   Map<GKInstance, String> entityToKey) throws Exception {
        if (entity.getSchemClass().isa(ReactomeJavaConstants.Complex))
            generateKeyForComplex(entity, entityToKey);
        else if (entity.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
            generateKeyForSet(entity, entityToKey);
        else
            generateKeyForAtomicEntity(entity, entityToKey);
    }
    
    private String generateKeyForComplex(GKInstance complex,
                                         Map<GKInstance, String> entityToString) throws Exception {
        String key = entityToString.get(complex);
        if (key != null)
            return key;
        StringBuilder builder = new StringBuilder();
        builder.append("Complex:");
        List<GKInstance> components = new ArrayList<GKInstance>();
        flatComplex(complex, components);
        if (components.size() == 0)
            builder.append(complex.getDisplayName());
        else {
            List<String> keys = new ArrayList<String>();
            for (GKInstance comp : components) {
                String key1 = null;
                if (comp.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
                    key1 = generateKeyForSet(comp, entityToString);
                else
                    key1 = generateKeyForAtomicEntity(comp, entityToString);
                keys.add(key1);
            }
            Collections.sort(keys);
            for (Iterator<String> it = keys.iterator(); it.hasNext();) {
                builder.append(it.next());
                if (it.hasNext())
                    builder.append(",");
            }
        }
        key = builder.toString();
        entityToString.put(complex, key);
        return key;
    }
    
    private void flatComplex(GKInstance complex, List<GKInstance> list) throws Exception {
        List comps = complex.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
        if (comps != null && comps.size() > 0) {
            for (Iterator it = comps.iterator(); it.hasNext();) {
                GKInstance comp = (GKInstance) it.next();
                SchemaClass cls = comp.getSchemClass();
                if (cls.isa(ReactomeJavaConstants.Complex))
                    flatComplex(comp, list);
                else
                    list.add(comp);
            }
        }
    }
    
    private void flatEntitySet(GKInstance set, List<GKInstance> list) throws Exception {
        List members = set.getAttributeValuesList(ReactomeJavaConstants.hasMember);
        if (members != null && members.size() > 0) {
            for (Iterator it = members.iterator(); it.hasNext();) {
                GKInstance member = (GKInstance) it.next();
                SchemaClass cls = member.getSchemClass();
                if (cls.isa(ReactomeJavaConstants.EntitySet))
                    flatEntitySet(member, list);
                else
                    list.add(member);
            }
        }
    }
    
    private String generateKeyForSet(GKInstance set,
                                     Map<GKInstance, String> entityToString) throws Exception {
        String key = entityToString.get(set);
        if (key != null)
            return key;
        StringBuilder builder = new StringBuilder();
        builder.append("Set:");
        List<GKInstance> members = new ArrayList<GKInstance>();
        flatEntitySet(set, members);
        if (members.size() == 0) {
            builder.append(set.getDisplayName());
        }
        else {
            List<String> keys = new ArrayList<String>();
            for (GKInstance member : members) {
                String key1 = generateKeyForAtomicEntity(member, entityToString);
                keys.add(key1);
            }
            Collections.sort(keys);
            for (Iterator<String> it = keys.iterator(); it.hasNext();) {
                builder.append(it.next());
                if (it.hasNext())
                    builder.append("|");
            }
        }
        key = builder.toString();
        entityToString.put(set, key);
        return key;
    }
        
    private String generateKeyForAtomicEntity(GKInstance entity,
                                             Map<GKInstance, String> entityToKey) throws Exception {
        String key = entityToKey.get(entity);
        if (key != null)
            return key;
        StringBuilder builder = new StringBuilder();
        boolean isDisplayNameUsed = false;
        // Use ids from ReferenceEntities
        if (entity.getSchemClass().isValidAttribute(ReactomeJavaConstants.referenceEntity)) {
            GKInstance refEntity = (GKInstance) entity.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            if (refEntity != null) {
                String identifier = (String) refEntity.getAttributeValue(ReactomeJavaConstants.identifier);
                if (identifier != null && identifier.length() > 0)
                    builder.append(identifier);
            }
        }
        // Otherwise use _displayName
        if (builder.length() == 0) {
            builder.append(entity.getDisplayName()); // Very weak. Should be avoided!
            isDisplayNameUsed = true;
        }
        // Have to consider modification
        if (entity.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasModifiedResidue)) {
            List modifications = entity.getAttributeValuesList(ReactomeJavaConstants.hasModifiedResidue);
            if (modifications != null && modifications.size() > 0) {
                for (Iterator it = modifications.iterator(); it.hasNext();) {
                    GKInstance modification = (GKInstance) it.next();
                    builder.append("[");
                    GKInstance chemical = (GKInstance) modification.getAttributeValue(ReactomeJavaConstants.modification);
                    if (chemical == null)
                        builder.append(modification.getDisplayName());
                    else
                        builder.append(chemical.getDisplayName());
                    Integer pos = (Integer) modification.getAttributeValue(ReactomeJavaConstants.coordinate);
                    if (pos != null)
                        builder.append(" at " ).append(pos);
                    builder.append("]");
                }
            }
        }
        // and compartment
        if (!isDisplayNameUsed && entity.getSchemClass().isValidAttribute(ReactomeJavaConstants.compartment)) {
            GKInstance compartment = (GKInstance) entity.getAttributeValue(ReactomeJavaConstants.compartment);
            if (compartment != null) 
                builder.append("[").append(compartment.getDisplayName()).append("]");
        }
        key = builder.toString();
        entityToKey.put(entity, 
                        key);
        return key;
    }
    
    @Test
    public void testMergeLocalProject() throws Exception {
        String fileName = "results/v2/TwoApoptosis.rtpj";
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        fileAdaptor.setSource(fileName);
        Collection reactions = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Reaction);
        List<GKInstance> before = new ArrayList<GKInstance>(reactions);
        System.out.println("Before: " + before.size());
        List<GKInstance> after = merge(before);
        System.out.println("After: " + after.size());
    }
    
    @Test
    @SuppressWarnings("unchecked")
    public void checkRedundancyOfReactions() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_plus_i_v2",
                                            "root",
                                            "macmysql01",
                                            3306);
        Collection reactions = prepareReactions(dba);
        List<GKInstance> set = new ArrayList<GKInstance>(reactions);
        System.out.println("Before merging: " + set.size());
        List<GKInstance> list = merge(set);
        System.out.println("After merging: " + list.size());
    }
    
    private Collection prepareReactions(MySQLAdaptor dba) throws Exception {
        // Load all reactions for analyzed
        Collection reactions = dba.fetchInstancesByClass(ReactomeJavaConstants.Reaction);
        Collection cas = dba.fetchInstancesByClass(ReactomeJavaConstants.CatalystActivity);
        Collection regulations = dba.fetchInstancesByClass(ReactomeJavaConstants.Regulation);
        Collection entities = dba.fetchInstancesByClass(ReactomeJavaConstants.EntityWithAccessionedSequence);
        // Load precedingEvent values
        SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Reaction);
        SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.input);
        dba.loadInstanceAttributeValues(reactions, att);
        att = cls.getAttribute(ReactomeJavaConstants.output);
        dba.loadInstanceAttributeValues(reactions, att);
        att = cls.getAttribute(ReactomeJavaConstants.catalystActivity);
        dba.loadInstanceAttributeValues(reactions, att);
        cls = dba.getSchema().getClassByName(ReactomeJavaConstants.CatalystActivity);
        att = cls.getAttribute(ReactomeJavaConstants.physicalEntity);
        dba.loadInstanceAttributeValues(cas, att);
        cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Regulation);
        att = cls.getAttribute(ReactomeJavaConstants.regulatedEntity);
        dba.loadInstanceAttributeValues(regulations, att);
        att = cls.getAttribute(ReactomeJavaConstants.regulator);
        dba.loadInstanceAttributeValues(regulations, att);
        cls = dba.getSchema().getClassByName(ReactomeJavaConstants.EntityWithAccessionedSequence);
        att = cls.getAttribute(ReactomeJavaConstants.referenceEntity);
        dba.loadInstanceAttributeValues(entities, att);
        Collection complexes = dba.fetchInstancesByClass(ReactomeJavaConstants.Complex);
        cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Complex);
        att = cls.getAttribute(ReactomeJavaConstants.hasComponent);
        dba.loadInstanceAttributeValues(complexes, att);
        return reactions;
    }
    
}
