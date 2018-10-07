package org.sbolstandard.AddGene2SBOL;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.net.URI;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.TimeZone;

import javax.xml.namespace.QName;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.sbolstandard.core2.AccessType;
import org.sbolstandard.core2.Annotation;
import org.sbolstandard.core2.Component;
import org.sbolstandard.core2.ComponentDefinition;
import org.sbolstandard.core2.GenericTopLevel;
import org.sbolstandard.core2.Implementation;
import org.sbolstandard.core2.SBOLConversionException;
import org.sbolstandard.core2.SBOLDocument;
import org.sbolstandard.core2.SBOLReader;
import org.sbolstandard.core2.SBOLValidate;
import org.sbolstandard.core2.SBOLValidationException;
import org.sbolstandard.core2.Sequence;
import org.sbolstandard.core2.SequenceAnnotation;
import org.sbolstandard.core2.SnapGene;
import org.sbolstandard.core2.TopLevel;
import org.synbiohub.frontend.SynBioHubException;
import org.synbiohub.frontend.SynBioHubFrontend;

public class AddGene2SBOL {
	
	private static String uriPrefix = "http://addgene.org/"; 
	private static String version = "1";
	private static String so = "http://identifiers.org/so/";
	private static String provNS = "http://www.w3.org/ns/prov#";
	private static String dcNS = "http://purl.org/dc/elements/1.1/";
	//private static String dcTermsNS = "http://purl.org/dc/terms/";
	private static String addgeneNS = "http://addgene.org/Terms/addgene#";
	private static String oboNS = "http://purl.obolibrary.org/obo/";

	static URI activityURI;
	static String createdDate;
	
	private static String createTimeString(long time1, long time2)
	{
		long minutes;
		long hours;
		long days;
		double secs = ((time2 - time1) / 1000000000.0);
		long seconds = ((time2 - time1) / 1000000000);
		secs = secs - seconds;
		minutes = seconds / 60;
		secs = seconds % 60 + secs;
		hours = minutes / 60;
		minutes = minutes % 60;
		days = hours / 24;
		hours = hours % 60;
		String time;
		String dayLabel;
		String hourLabel;
		String minuteLabel;
		String secondLabel;
		if (days == 1)
		{
			dayLabel = " day ";
		}
		else
		{
			dayLabel = " days ";
		}
		if (hours == 1)
		{
			hourLabel = " hour ";
		}
		else
		{
			hourLabel = " hours ";
		}
		if (minutes == 1)
		{
			minuteLabel = " minute ";
		}
		else
		{
			minuteLabel = " minutes ";
		}
		if (seconds == 1)
		{
			secondLabel = " second";
		}
		else
		{
			secondLabel = " seconds";
		}
		if (days != 0)
		{
			time = days + dayLabel + hours + hourLabel + minutes + minuteLabel + secs + secondLabel;
		}
		else if (hours != 0)
		{
			time = hours + hourLabel + minutes + minuteLabel + secs + secondLabel;
		}
		else if (minutes != 0)
		{
			time = minutes + minuteLabel + secs + secondLabel;
		}
		else
		{
			time = secs + secondLabel;
		}
		return time;
	}
	
	private static void addgeneStringAnnotation(List<Annotation> annotations, JSONObject object, String tag) throws SBOLValidationException {
		String annotationValue = object.get(tag)!=null?object.get(tag).toString():null;
		if (annotationValue!=null && !annotationValue.equals("")) {
			Annotation annotation = new Annotation(new QName(addgeneNS,tag,"ag"), annotationValue);
			annotations.add(annotation);
		}
	}

	private static void addgeneStringAnnotationArray(List<Annotation> annotations, JSONObject object, String tag) throws SBOLValidationException {
		JSONArray annotationValues = (JSONArray)object.get(tag);
		for (Object annotationValue : annotationValues) {
			Annotation annotation = new Annotation(new QName(addgeneNS,tag,"ag"), (String)annotationValue);
			annotations.add(annotation);
		}
	}
	
	private static void createAddgeneStringAnnotation(TopLevel topLevel,JSONObject plasmid, String tag) throws SBOLValidationException {
		String annotation = plasmid.get(tag)!=null?plasmid.get(tag).toString():null;
		if (annotation!=null && !annotation.equals("")) {
			topLevel.createAnnotation(new QName(addgeneNS,tag,"ag"), annotation);
		}
	}
		
	private static void createAddgeneStringAnnotation(SBOLDocument document, Implementation implementation, JSONObject plasmid, String tag) 
			throws SBOLValidationException {
		createAddgeneStringAnnotation(implementation,plasmid,tag);
		if (implementation.isSetBuilt()) {
			createArticleAnnotations((TopLevel)implementation.getBuilt(),plasmid);
			ComponentDefinition cd = (ComponentDefinition)implementation.getBuilt();
			createAddgeneStringAnnotation(cd,plasmid,tag);
			Sequence seq = cd.getSequenceByEncoding(Sequence.IUPAC_DNA);
			if (seq!=null) {
				createAddgeneStringAnnotation(seq,plasmid,tag);
			}
		} 
		for (URI wasDerivedFrom : implementation.getWasDerivedFroms()) {
			if (!wasDerivedFrom.equals(implementation.getBuiltURI())) {
				ComponentDefinition cd = document.getComponentDefinition(wasDerivedFrom);
				createAddgeneStringAnnotation(cd,plasmid,tag);
				Sequence seq = cd.getSequenceByEncoding(Sequence.IUPAC_DNA);
				if (seq!=null) {
					createAddgeneStringAnnotation(seq,plasmid,tag);
				}
			}
		}
	}
		
	private static void createAddgeneStringAnnotationArray(TopLevel topLevel,JSONObject object, String tag) throws SBOLValidationException {
		JSONArray annotationValues = (JSONArray)object.get(tag);
		for (Object annotationValue : annotationValues) {
			topLevel.createAnnotation(new QName(addgeneNS,tag,"ag"), (String)annotationValue);
		}
	}
	
	private static void createSequence(SBOLDocument document, ComponentDefinition cd,
			JSONObject sequences, String tag, String idSuffix) throws SBOLValidationException {
		JSONArray seqObj = (JSONArray)sequences.get(tag);
		int j = 0;
		for (Object seq : seqObj) {
			Sequence s = document.createSequence(cd.getDisplayId()+idSuffix+j, version, ((String)seq).trim(), Sequence.IUPAC_DNA);
			cd.addSequence(s);
			j++;
		}
	}
	
	private static void createArticleAnnotations(TopLevel topLevel, JSONObject plasmid) throws SBOLValidationException {
		JSONObject article = (JSONObject)plasmid.get("article");
		if (article != null) {
			if (article.get("pubmed_id")!=null) {
				String pubmedId = article.get("pubmed_id").toString();
				topLevel.createAnnotation(new QName(oboNS,"OBI_0001617","obo"), pubmedId);
			}
			if (article.get("doi")!=null) {
				String doi = article.get("doi").toString();
				topLevel.createAnnotation(new QName(oboNS,"OBI_0002110","obo"), doi);
			}
			if (article.get("id")!=null) {
				String id = article.get("id").toString();
				topLevel.createAnnotation(new QName(addgeneNS,"articleId","ag"), id);
			}
			if (article.get("url")!=null) {
				URI articleUrl = URI.create(article.get("url").toString());
				topLevel.createAnnotation(new QName(addgeneNS,"articleUrl","ag"), articleUrl);
			}
		}
	}
	
	private static void createArticleAnnotations(SBOLDocument document, Implementation implementation, JSONObject plasmid) throws SBOLValidationException {
		createArticleAnnotations((TopLevel)implementation,plasmid);
		if (implementation.isSetBuilt()) {
			createArticleAnnotations((TopLevel)implementation.getBuilt(),plasmid);
			ComponentDefinition cd = (ComponentDefinition)implementation.getBuilt();
			createArticleAnnotations((TopLevel)cd,plasmid);
			Sequence seq = cd.getSequenceByEncoding(Sequence.IUPAC_DNA);
			if (seq!=null) {
				createArticleAnnotations((TopLevel)seq,plasmid);
			}
		} 
		for (URI wasDerivedFrom : implementation.getWasDerivedFroms()) {
			if (!wasDerivedFrom.equals(implementation.getBuiltURI())) {
				ComponentDefinition cd = document.getComponentDefinition(wasDerivedFrom);
				createArticleAnnotations((TopLevel)cd,plasmid);
				Sequence seq = cd.getSequenceByEncoding(Sequence.IUPAC_DNA);
				if (seq!=null) {
					createArticleAnnotations((TopLevel)seq,plasmid);
				}
			}
		}
	}
	
	private static void createCreatorAnnotations(TopLevel topLevel, JSONObject plasmid) throws SBOLValidationException {
		JSONArray pis = (JSONArray)plasmid.get("pi");
		for (Object pi : pis) {
			topLevel.createAnnotation(new QName(dcNS,"creator","dc"), (String)pi);
		}
	}
		
	private static void createCreatorAnnotations(SBOLDocument document, Implementation implementation, JSONObject plasmid) throws SBOLValidationException {
		createCreatorAnnotations((TopLevel)implementation,plasmid);
		if (implementation.isSetBuilt()) {
			createCreatorAnnotations((TopLevel)implementation.getBuilt(),plasmid);
			ComponentDefinition cd = (ComponentDefinition)implementation.getBuilt();
			createCreatorAnnotations((TopLevel)cd,plasmid);
			Sequence seq = cd.getSequenceByEncoding(Sequence.IUPAC_DNA);
			if (seq!=null) {
				createCreatorAnnotations((TopLevel)seq,plasmid);
			}
		} 
		for (URI wasDerivedFrom : implementation.getWasDerivedFroms()) {
			if (!wasDerivedFrom.equals(implementation.getBuiltURI())) {
				ComponentDefinition cd = document.getComponentDefinition(wasDerivedFrom);
				createCreatorAnnotations((TopLevel)cd,plasmid);
				Sequence seq = cd.getSequenceByEncoding(Sequence.IUPAC_DNA);
				if (seq!=null) {
					createCreatorAnnotations((TopLevel)seq,plasmid);
				}
			}
		}
	}
	
	private static void createCloningAnnotations(TopLevel topLevel, JSONObject plasmid) throws SBOLValidationException {
		JSONObject cloning = (JSONObject)plasmid.get("cloning");
		List<Annotation> annotations = new ArrayList<Annotation>();
		addgeneStringAnnotation(annotations,cloning,"backbone");
		addgeneStringAnnotation(annotations,cloning,"backbone_mutation");
		addgeneStringAnnotation(annotations,cloning,"backbone_origin");
		addgeneStringAnnotation(annotations,cloning,"backbone_size");
		addgeneStringAnnotation(annotations,cloning,"promoter");
		addgeneStringAnnotation(annotations,cloning,"sequencing_primer_3");
		addgeneStringAnnotation(annotations,cloning,"sequencing_primer_5");
		addgeneStringAnnotationArray(annotations,cloning,"vector_types");
		topLevel.createAnnotation(new QName(addgeneNS,"cloning","ag"), new QName(addgeneNS,"Cloning","ag"), 
				URI.create(topLevel.getPersistentIdentity().toString()+"/cloning"), annotations);
	}
	
	private static void createTagAnnotation(List<Annotation> tagAnnotations, JSONObject object, 
			String nestedUriPrefix) throws SBOLValidationException {
		JSONArray tags = (JSONArray)object.get("tags");
		int i = 0;
		for (Object tag : tags) {
			List<Annotation> annotations = new ArrayList<Annotation>();
			addgeneStringAnnotation(annotations,(JSONObject)tag,"location");
			addgeneStringAnnotation(annotations,(JSONObject)tag,"tag");
			Annotation annotation = new Annotation(new QName(addgeneNS,"tagElement","ag"), new QName(addgeneNS,"TagElement","ag"), 
					URI.create(nestedUriPrefix+"/tagElement"+i), annotations);
			i++;
			tagAnnotations.add(annotation);
		}	
	}
	
	private static void createInsertCloningAnnotations(List<Annotation> insertAnnotations, 
			JSONObject insert, int i, String nestedUriPrefix) throws SBOLValidationException {
		JSONObject cloning = (JSONObject)insert.get("cloning");
		if (cloning==null) return;
		List<Annotation> annotations = new ArrayList<Annotation>();
		addgeneStringAnnotation(annotations,(JSONObject)cloning,"clone_method");
		addgeneStringAnnotation(annotations,(JSONObject)cloning,"cloning_site_3");
		addgeneStringAnnotation(annotations,(JSONObject)cloning,"cloning_site_5");
		addgeneStringAnnotation(annotations,(JSONObject)cloning,"promoter");
		addgeneStringAnnotation(annotations,(JSONObject)cloning,"sequencing_primer_3");
		addgeneStringAnnotation(annotations,(JSONObject)cloning,"sequencing_primer_5");
		addgeneStringAnnotation(annotations,(JSONObject)cloning,"site_3_destroyed");
		addgeneStringAnnotation(annotations,(JSONObject)cloning,"site_5_destroyed");
		Annotation annotation = new Annotation(new QName(addgeneNS,"cloning","ag"), new QName(addgeneNS,"Cloning","ag"), 
				URI.create(nestedUriPrefix+"/insert"+i+"/cloning"), annotations);
		i++;
		insertAnnotations.add(annotation);
	}
	
	private static void createInsertEntrezGeneAnnotations(List<Annotation> insertAnnotations, 
			JSONObject insert, int insertNum, String nestedUriPrefix) throws SBOLValidationException {
		JSONArray entrezGenes = (JSONArray)insert.get("entrez_gene");
		int i = 0;
		for (Object entrezGene : entrezGenes) {
			List<Annotation> annotations = new ArrayList<Annotation>();
			addgeneStringAnnotation(annotations,(JSONObject)entrezGene,"aliases");
			addgeneStringAnnotation(annotations,(JSONObject)entrezGene,"gene");
			addgeneStringAnnotation(annotations,(JSONObject)entrezGene,"id");
			Annotation annotation = new Annotation(new QName(addgeneNS,"entrezGene","ag"), new QName(addgeneNS,"EntrezGene","ag"), 
					URI.create(nestedUriPrefix+"/insert"+ insertNum +"/entrezGene"+i), annotations);
			i++;
			insertAnnotations.add(annotation);
		}	
	}
	
	private static void createInsertSpeciesAnnotations(List<Annotation> insertAnnotations, JSONObject insert) throws SBOLValidationException {
		JSONArray speciesArrayOuter = (JSONArray)insert.get("species");
		for (Object speciesArrayInner : speciesArrayOuter) {
			for (Object species : (JSONArray)speciesArrayInner) {
				if (species!=null) {
					Annotation annotation = new Annotation(new QName(addgeneNS,"species","ag"), species.toString());
					insertAnnotations.add(annotation);
				}
			}
		}
	}
	
	private static void createInsertAnnotations(TopLevel topLevel, JSONObject plasmid) throws SBOLValidationException {
		JSONArray inserts = (JSONArray)plasmid.get("inserts");
		int i = 0;
		for (Object insert : inserts) {
			List<Annotation> annotations = new ArrayList<Annotation>();
			addgeneStringAnnotationArray(annotations,(JSONObject)insert,"alt_names");
			createInsertCloningAnnotations(annotations,(JSONObject)insert,i,topLevel.getPersistentIdentity().toString());
			createInsertEntrezGeneAnnotations(annotations,(JSONObject)insert,i,topLevel.getPersistentIdentity().toString());
			createInsertSpeciesAnnotations(annotations,(JSONObject)insert);
			addgeneStringAnnotationArray(annotations,(JSONObject)insert,"genbank_ids");
			addgeneStringAnnotation(annotations,(JSONObject)insert,"mutation");
			addgeneStringAnnotation(annotations,(JSONObject)insert,"name");
			addgeneStringAnnotation(annotations,(JSONObject)insert,"shRNA_sequence");
			addgeneStringAnnotation(annotations,(JSONObject)insert,"size");
			createTagAnnotation(annotations,(JSONObject)insert,topLevel.getPersistentIdentity().toString()+"/insert"+i);
			topLevel.createAnnotation(new QName(addgeneNS,"insert","ag"), new QName(addgeneNS,"Insert","ag"), 
					URI.create(topLevel.getPersistentIdentity().toString()+"/insert"+i), annotations);
			i++;
		}
	}
	
	private static void detectFeatures(SBOLDocument document, String displayId, String name, JSONObject plasmid) throws SBOLValidationException {
		JSONObject sequences = (JSONObject) plasmid.get("sequences");		
		JSONArray seqObj = (JSONArray)sequences.get("public_addgene_full_sequences");
		if (seqObj.size()==0) {
			seqObj = (JSONArray)sequences.get("public_user_full_sequences");
			if (seqObj.size()==0) {
				seqObj = (JSONArray)sequences.get("public_addgene_partial_sequences");
				if (seqObj.size()==0) {
					seqObj = (JSONArray)sequences.get("public_user_partial_sequences");
				}
			}
		}
		for (Object seq : seqObj) {
			SBOLDocument doc = SnapGene.detectFeatures(((String)seq).trim(),uriPrefix,displayId,version,
					"/tmp/addGene/"+displayId+".png",name);
			if (doc!=null) {
				document.createCopy(doc);
			}
			return;
		}
	}
	
	private static void statistics(JSONArray plasmids) {
		
		int au_equal = 0;
		int au_notEqual = 0;
		int au_total = 0;
		int a_equal = 0;
		int a_notEqual = 0;
		int a_total = 0;
		int ap_found = 0;
		int ap_notFound = 0;
		int ap_total = 0;
		int u_equal = 0;
		int u_notEqual = 0;
		int u_total = 0;
		int up_found = 0;
		int up_notFound = 0;
		int up_total = 0;
		
		for (Object p : plasmids) {
			JSONObject aPlasmid = (JSONObject) p;
			String aDisplayId = "addgene"+aPlasmid.get("id");
			JSONObject allSequences = (JSONObject) aPlasmid.get("sequences");
			JSONArray addSeqObj = (JSONArray)allSequences.get("public_addgene_full_sequences");
			JSONArray userSeqObj = (JSONArray)allSequences.get("public_user_full_sequences");
			JSONArray addPSeqObj = (JSONArray)allSequences.get("public_addgene_partial_sequences");
			JSONArray userPSeqObj = (JSONArray)allSequences.get("public_user_partial_sequences");
			if (userSeqObj.size()>0 && addSeqObj.size() > 0) {
				String seqString = ((String)addSeqObj.get(0)).trim();
				System.out.println(aDisplayId + " : number of addgene_full_seqenences = " + addSeqObj.size());
				System.out.println(aDisplayId + " : number of user_full_seqenences = " + userSeqObj.size());
				au_total++;
				for (int j = 0; j < userSeqObj.size(); j++) {
					String seqString2 = ((String)userSeqObj.get(j)).trim();
					if (seqString.equals(seqString2)) {
						System.out.println("  "+j+" : equal to 0");
						au_equal++;
					} else {
						System.out.println("  "+j+" : not equal to 0");
						au_notEqual++;
					}
				}
			} 
			if (userSeqObj.size()>1) {
				String seqString = ((String)userSeqObj.get(0)).trim();
				System.out.println(aDisplayId + " : number of user_full_seqenences = " + userSeqObj.size());
				u_total++;
				for (int j = 1; j < userSeqObj.size(); j++) {
					String seqString2 = ((String)userSeqObj.get(j)).trim();
					if (seqString.equals(seqString2)) {
						System.out.println("  "+j+" : equal to 0");
						u_equal++;
					} else {
						System.out.println("  "+j+" : not equal to 0");
						u_notEqual++;
					}
				}
			} 
			if (addSeqObj.size() > 1) {
				String seqString = ((String)addSeqObj.get(0)).trim();
				System.out.println(aDisplayId + " : number of addgene_full_seqenences = " + addSeqObj.size());
				a_total++;
				for (int j = 1; j < addSeqObj.size(); j++) {
					String seqString2 = ((String)addSeqObj.get(j)).trim();
					if (seqString.equals(seqString2)) {
						System.out.println("  "+j+" : equal to 0");
						a_equal++;
					} else {
						System.out.println("  "+j+" : not equal to 0");
						a_notEqual++;
					}
				}
			} 
			if (addSeqObj.size() > 0 && addPSeqObj.size() > 0) {
				System.out.println(aDisplayId + " : number of addgene_partial_seqenences = " + addPSeqObj.size());
				ap_total++;
				String seqString = ((String)addSeqObj.get(0)).trim();
				for (int j = 0; j < addPSeqObj.size(); j++) {
					String seqString2 = ((String)addPSeqObj.get(j)).trim();
					if (seqString.contains(seqString2)) {
						System.out.println("  "+j+" : found partial");
						ap_found++;
					} else {
						System.out.println("  "+j+" : not found partial");
						ap_notFound++;
					}
				}				
			}
			if (userSeqObj.size() > 0 && userPSeqObj.size() > 0) {
				System.out.println(aDisplayId + " : number of user_partial_seqenences = " + userPSeqObj.size());
				up_total++;
				String seqString = ((String)userSeqObj.get(0)).trim();
				for (int j = 0; j < userPSeqObj.size(); j++) {
					String seqString2 = ((String)userPSeqObj.get(j)).trim();
					if (seqString.contains(seqString2)) {
						System.out.println("  "+j+" : found partial");
						up_found++;
					} else {
						System.out.println("  "+j+" : not found partial");
						up_notFound++;
					}
				}				
			}
			System.out.println("Total : (addgene_full) " + a_total + " with duplicates " + 
					a_equal + " equal duplicates " + a_notEqual + " not equal duplicates");
			System.out.println("Total : (user_full) " + u_total + " with duplicates " + 
					u_equal + " equal duplicates " + u_notEqual + " not equal duplicates");
			System.out.println("Total : (both addgene_full and user_full) " + au_total + " with both " + 
					au_equal + " equal " + au_notEqual + " not equal");
			System.out.println("Total : (both addgene_full and addgene_partial) " + ap_total + " with both " + 
					(ap_found + ap_notFound) + " total partials " + ap_found + " found " + ap_notFound + " not found");
			System.out.println("Total : (both user_full and user_partial) " + up_total + " with both " + 
					(up_found + up_notFound) + " total partials " + up_found + " found " + up_notFound + " not found");
		}
	}
	
	private static void createBenData() {
		// Conversion of Ben's data, belongs in different project
//		SynBioHubFrontend sbh = new SynBioHubFrontend("https://synbiohub.utah.edu");
//		sbh.login("myers@ece.utah.edu", "MaWen69!");
//		sbh.createCollection("BenSplitRecomb", "1", "Ben's Split Recombinases", "Ben's Split Recombinases", "", true);
//		File folder = new File("/Users/myers/Downloads/BenData/Split Recombinases for Addgene/");
//		File[] listOfFiles = folder.listFiles();
//		for (File file : listOfFiles) {
//			try {
//				System.out.println(file.getName());
//				SBOLDocument finalDoc = new SBOLDocument();
//				String uriPrefix = "http://ben.data.org/";
//				SBOLReader.setURIPrefix(uriPrefix);
//				SBOLReader.setDisplayId("");
//				SBOLDocument doc = SBOLReader.read(file);
//				ComponentDefinition cd = doc.getRootComponentDefinitions().iterator().next();
//				Sequence seq = cd.getSequenceByEncoding(Sequence.IUPAC_DNA);
//				finalDoc = SynBioHubFrontend.detectFeatures("http://tang.ece.utah.edu/examples/pages/acceptNewText.php",
//						"http://tang.ece.utah.edu/dnafiles/",seq.getElements(),uriPrefix,cd.getDisplayId(),cd.getVersion());
//				sbh.addToCollection(URI.create("https://synbiohub.utah.edu/user/myers/BenSplitRecomb/BenSplitRecomb_collection/1"), false, finalDoc);
//			} catch (Exception e) {
//				e.printStackTrace();
//			}
//		}
//		if (true) return;
	}
	
	private static ComponentDefinition createComponentDefinition(SBOLDocument document, String sequence, String displayId, 
			String name,boolean circular) throws SBOLValidationException {
		if (sequence != null) {
			SBOLDocument doc = SnapGene.detectFeatures(sequence,uriPrefix,displayId,version,"/tmp/addGene/"+displayId+".png",name,circular);
			if (doc!=null) {
				document.createCopy(doc);
			}
		}
		ComponentDefinition cd = document.getComponentDefinition(displayId, version);
		if (cd==null) {
			cd = document.createComponentDefinition(displayId, version, ComponentDefinition.DNA);
			if (sequence != null) {
				document.createSequence(displayId+"_seq", version, sequence, Sequence.IUPAC_DNA);
			}
		}
		return cd;
	}
	
	private static boolean sequenceMismatch(JSONArray seqObj) {
		String seqString = ((String)seqObj.get(0)).trim();
		for (int j = 1; j < seqObj.size(); j++) {
			String seqString2 = ((String)seqObj.get(j)).trim();
			if (!seqString.equals(seqString2)) {
				return true;
			}
		}
		return false;
	}
	
	public static void main( String[] args ) throws FileNotFoundException, IOException, ParseException, SynBioHubException, SBOLConversionException, SBOLValidationException 
    {
		// Create an SBOLDocument
		SBOLDocument document = new SBOLDocument(); 
		document.setDefaultURIprefix(uriPrefix); 
		document.setComplete(true); 
		document.setCreateDefaults(true);
		//URI activity = null;
		SynBioHubFrontend sbh = new SynBioHubFrontend(args[4],args[3]);
		// Create collection
		//System.out.println(args[0]); // login
		//System.out.println(args[1]); // password
		//System.out.println(args[2]); // JSON file
		//System.out.println(args[3]); // URI prefix
		//System.out.println(args[4]); // URL
		try {
			sbh.login(args[0], args[1]);
		}
		catch (SynBioHubException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
			return;
		}
		int start = 0;
		
		try {
			// Create an Activity
			GenericTopLevel genericTopLevel = document.createGenericTopLevel("addgene2sbol", version, 
					new QName(provNS, "Activity", "prov"));
			genericTopLevel.getIdentity();
			genericTopLevel.setName("Addgene to SBOL conversion");
			genericTopLevel.setDescription("Conversion of the Addgene plasmids and metadata to SBOL2");
			genericTopLevel.createAnnotation(new QName(dcNS,"creator","dc"), "Chris J. Myers");
			TimeZone tz = TimeZone.getTimeZone("UTC");
			DateFormat df = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss.SSS'Z'");
			df.setTimeZone(tz);
			createdDate = df.format(new Date());
			genericTopLevel.createAnnotation(new QName(provNS,"endedAtTime","prov"), createdDate);
			activityURI = genericTopLevel.getIdentity();

			if (start==0)
				sbh.createCollection("AddgenePlasmids", "1", "Addgene Plasmids", "These are the Addgene plasmids", "", true, document);
			//document.write(System.out);
		} catch (Exception e) {
			e.printStackTrace();
		}
		// Parse JSON
		JSONParser parser = new JSONParser();
		JSONObject o = (JSONObject) parser.parse(new FileReader(args[2]));

		JSONArray plasmids = (JSONArray) o.get("plasmids");
		int i = 0;
		int size = plasmids.size();
		int success = 0;
		int failure = 0;
		int skip = 0;
		int a_mismatch = 0;
		int u_mismatch = 0;
		long time0, time1, time2;
		time0 = System.nanoTime();
		
		for (Object p : plasmids) {
			time1 = System.nanoTime();
			i++;
			if (i < start) continue;
			document = new SBOLDocument(); 
			document.setDefaultURIprefix(uriPrefix); 
			//document.setComplete(true); 
			document.setCreateDefaults(true);
			JSONObject plasmid = (JSONObject) p;
			String displayId = "addgene"+plasmid.get("id");
			String name = plasmid.get("name").toString();
			try {
				JSONObject sequences = (JSONObject) plasmid.get("sequences");
				Implementation imp = document.getImplementation(displayId+"_plasmid", version);
				ComponentDefinition cd = null;
				ComponentDefinition ucd = null;
				if (imp==null) {
					imp = document.createImplementation(displayId+"_plasmid", version);
					JSONArray addFullSeqObj = (JSONArray)sequences.get("public_addgene_full_sequences");
					JSONArray userFullSeqObj = (JSONArray)sequences.get("public_user_full_sequences");
					if (addFullSeqObj.size()==1 && userFullSeqObj.size()==0) {
						cd = createComponentDefinition(document,((String)addFullSeqObj.get(0)).trim(),displayId,name,true);
						imp.setBuilt(cd.getIdentity());
						imp.addWasDerivedFrom(cd.getIdentity());
					} else if (addFullSeqObj.size()==0 && userFullSeqObj.size()==1) {
						cd = createComponentDefinition(document,((String)userFullSeqObj.get(0)).trim(),displayId,name,true);
						imp.setBuilt(cd.getIdentity());
						imp.addWasDerivedFrom(cd.getIdentity());
					} else if (addFullSeqObj.size()==1 && userFullSeqObj.size()==1 &&
							((String)addFullSeqObj.get(0)).trim().equals((String)userFullSeqObj.get(0))) {
						cd = createComponentDefinition(document,((String)addFullSeqObj.get(0)).trim(),displayId,name,true);
						imp.setBuilt(cd.getIdentity());
						imp.addWasDerivedFrom(cd.getIdentity());
					} else if (addFullSeqObj.size()==1 && userFullSeqObj.size()==1 &&
							!((String)addFullSeqObj.get(0)).trim().equals((String)userFullSeqObj.get(0))) {
						cd = createComponentDefinition(document,((String)addFullSeqObj.get(0)).trim(),displayId,name,true);
						imp.setBuilt(cd.getIdentity());
						ucd = createComponentDefinition(document,((String)userFullSeqObj.get(0)).trim(),displayId+"_user",name,true);
						imp.addWasDerivedFrom(ucd.getIdentity());
					} else if (addFullSeqObj.size() > 1) {
						if (sequenceMismatch(addFullSeqObj)) {
							a_mismatch++;
							time2 = System.nanoTime();
							String time = createTimeString(time1, time2);
							String totalTime = createTimeString(time0, time2);
							System.out.println(i + " out of " + size + ":"+displayId+" MULTIPLE ADDGENE FULL MISMATCH "+ a_mismatch + " " +time + " (" + totalTime +")");
							//System.err.println(i + " out of " + size + ":"+displayId+" SKIP "+ a_mismatch + " " +time + " (" + totalTime +")");
							continue;
						} else {
							cd = createComponentDefinition(document,((String)addFullSeqObj.get(0)).trim(),displayId,name,true);
							imp.setBuilt(cd.getIdentity());
							imp.addWasDerivedFrom(cd.getIdentity());
						}
					} else if (userFullSeqObj.size() > 1) {
						if (sequenceMismatch(userFullSeqObj)) {
							u_mismatch++;
							time2 = System.nanoTime();
							String time = createTimeString(time1, time2);
							String totalTime = createTimeString(time0, time2);
							System.out.println(i + " out of " + size + ":"+displayId+" MULTIPLE USER FULL MISMATCH "+ u_mismatch + " " +time + " (" + totalTime +")");
							//System.err.println(i + " out of " + size + ":"+displayId+" MULTIPLE USER FULL MISMATCH "+ u_mismatch + " " +time + " (" + totalTime +")");
							continue;
						} else {
							cd = createComponentDefinition(document,((String)userFullSeqObj.get(0)).trim(),displayId,name,true);
							imp.setBuilt(cd.getIdentity());
							imp.addWasDerivedFrom(cd.getIdentity());
						}
					} else if (addFullSeqObj.size()==0 && userFullSeqObj.size()==0) {
						cd = createComponentDefinition(document,null,displayId,name,true);
						imp.setBuilt(cd.getIdentity());
						imp.addWasDerivedFrom(cd.getIdentity());
					} else {
						skip++;
						time2 = System.nanoTime();
						String time = createTimeString(time1, time2);
						String totalTime = createTimeString(time0, time2);
			        	System.out.println(i + " out of " + size + ":"+displayId+" SKIP "+ skip + " " +time + " (" + totalTime +")");
			        	//System.err.println(i + " out of " + size + ":"+displayId+" SKIP "+ skip + " " +time + " (" + totalTime +")");
						continue;
					}
				}
				// TODO: partials
				imp.addWasGeneratedBy(URI.create(args[3]+"/user/myers/AddgenePlasmids/addgene2sbol/1"));
				cd.addWasGeneratedBy(URI.create(args[3]+"/user/myers/AddgenePlasmids/addgene2sbol/1"));
				if (ucd != null) {
					ucd.addWasGeneratedBy(URI.create(args[3]+"/user/myers/AddgenePlasmids/addgene2sbol/1"));
				}
				createArticleAnnotations(document,imp,plasmid);
				createAddgeneStringAnnotation(imp,plasmid,"bacterial_resistance");
				createCloningAnnotations(imp,plasmid);
				createAddgeneStringAnnotation(imp,plasmid,"growth_notes");
				createAddgeneStringAnnotation(imp,plasmid,"growth_strain");
				createAddgeneStringAnnotation(imp,plasmid,"growth_temp");
				createInsertAnnotations(imp,plasmid);
				createAddgeneStringAnnotation(imp,plasmid,"origin");
				createCreatorAnnotations(document,imp,plasmid);
				createAddgeneStringAnnotation(imp,plasmid,"plasmid_copy");
				createAddgeneStringAnnotationArray(imp,plasmid,"resistance_markers");
				List<Annotation> tagAnnotations = new ArrayList<Annotation>();
				createTagAnnotation(tagAnnotations,plasmid,imp.getPersistentIdentity().toString());
				if (tagAnnotations.size() > 0) {
					imp.createAnnotation(new QName(addgeneNS,"tags","ag"), new QName(addgeneNS,"Tags","ag"), 
							URI.create(imp.getPersistentIdentity()+"/tag"), tagAnnotations);
				}
				createAddgeneStringAnnotation(document,imp,plasmid,"terms");
				imp.addWasDerivedFrom(URI.create(plasmid.get("url").toString()));
				imp.setName(name);
				cd.setName(name);
				if (ucd!=null) {
					ucd.setName(name);
				}
				cd.addRole(URI.create(so + "SO:0000155"));
				JSONArray addPartialSeqObj = (JSONArray)sequences.get("public_addgene_partial_sequences");
				for (int j = 0; j < addPartialSeqObj.size(); j++) {
					String seqString = ((String)addPartialSeqObj.get(j)).trim();
					ComponentDefinition subCD = createComponentDefinition(document,seqString,displayId+"_addgene_partial"+j,null,false);
					Component comp = cd.createComponent(displayId+"_addgene_partial"+j, AccessType.PRIVATE, subCD.getIdentity());
					Sequence fullSeq = cd.getSequenceByEncoding(Sequence.IUPAC_DNA);
					if (fullSeq!=null) {
						String fullSeqStr = fullSeq.getElements();
						if (fullSeqStr.indexOf(seqString) >= 0) {
							SequenceAnnotation sa = cd.createSequenceAnnotation(displayId+"_addgene_partial_annotation"+j, "range", fullSeqStr.indexOf(seqString)+1, seqString.length());
							sa.setComponent(comp.getIdentity());
						}
					}
				}			
				JSONArray userPartialSeqObj = (JSONArray)sequences.get("public_user_partial_sequences");
				for (int j = 0; j < userPartialSeqObj.size(); j++) {
					String seqString = ((String)userPartialSeqObj.get(j)).trim();
					ComponentDefinition subCD = createComponentDefinition(document,seqString,displayId+"_user_partial"+j,null,false);
					ComponentDefinition topCD = cd;
					if (ucd!=null) {
						topCD = ucd;
					}
					Component comp = topCD.createComponent(displayId+"_user_partial"+j, AccessType.PRIVATE, subCD.getIdentity());
					Sequence fullSeq = topCD.getSequenceByEncoding(Sequence.IUPAC_DNA);
					if (fullSeq!=null) {
						String fullSeqStr = fullSeq.getElements();
						if (fullSeqStr.indexOf(seqString) >= 0) {
							SequenceAnnotation sa = topCD.createSequenceAnnotation(displayId+"_user_partial_annotation"+j, "range", fullSeqStr.indexOf(seqString)+1, seqString.length());
							sa.setComponent(comp.getIdentity());
						}
					}
				}
				
				SBOLValidate.validateSBOL(document,false,true,false);
				if (SBOLValidate.getNumErrors()>0) {
					failure++;
					time2 = System.nanoTime();
					String time = createTimeString(time1, time2);
					String totalTime = createTimeString(time0, time2);
					System.out.println(i + " out of " + size + ":"+displayId+" FAILURE "+ failure + " " +time + " (" + totalTime +")");
					System.err.println(i + " out of " + size + ":"+displayId+" FAILURE "+ failure + " " +time + " (" + totalTime +")");
					for (String error : SBOLValidate.getErrors()) {
						System.err.println(error);
					}
				} else {   
					// Upload to SynBioHub
					sbh.addToCollection(URI.create(args[3]+"/user/myers/AddgenePlasmids/AddgenePlasmids_collection/1"), true, document);
					File pngFile = new File("/tmp/addGene/"+cd.getDisplayId()+".png");
					if (pngFile.exists()) {
						sbh.attachFile(URI.create(args[3]+"/user/myers/AddgenePlasmids/"+cd.getDisplayId()+"/"+version), pngFile);
					}
					success++;
					time2 = System.nanoTime();
					String time = createTimeString(time1, time2);
					String totalTime = createTimeString(time0, time2);
		        	System.out.println(i + " out of " + size + ":"+displayId+" SUCCESS "+ success + " " +time + " (" + totalTime +")");
					//document.write(System.out);
					//System.exit(0);
				}
			} catch (Exception e) {
				failure++;
				time2 = System.nanoTime();
				String time = createTimeString(time1, time2);
				String totalTime = createTimeString(time0, time2);
	        	System.out.println(i + " out of " + size + ":"+displayId+" FAILURE "+ failure + " " +time + " (" + totalTime +")");
	        	System.err.println(i + " out of " + size + ":"+displayId+" FAILURE "+ failure + " " +time + " (" + totalTime +")");
	        	e.printStackTrace(System.err);
			}
		}
		//document.write(System.out);
        
        // Validate
    }
}
