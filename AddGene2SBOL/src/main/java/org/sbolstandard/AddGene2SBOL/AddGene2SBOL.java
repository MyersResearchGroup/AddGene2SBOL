package org.sbolstandard.AddGene2SBOL;

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

//import org.joda.time.DateTime;
//import org.joda.time.format.DateTimeFormatter;
//import org.joda.time.format.ISODateTimeFormat;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.sbolstandard.core2.Annotation;
import org.sbolstandard.core2.ComponentDefinition;
import org.sbolstandard.core2.GenericTopLevel;
import org.sbolstandard.core2.SBOLDocument;
import org.sbolstandard.core2.SBOLValidate;
import org.sbolstandard.core2.SBOLValidationException;
import org.sbolstandard.core2.Sequence;
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
	private static String addGeneNS = "http://addgene.org/Terms/addGene#";
	private static String oboNS = "http://purl.obolibrary.org/obo/";

	static URI activityURI;
	static String createdDate;
	
	private static void addGeneStringAnnotation(List<Annotation> annotations, JSONObject object, String tag) throws SBOLValidationException {
		String annotationValue = object.get(tag)!=null?object.get(tag).toString():null;
		if (annotationValue!=null && !annotationValue.equals("")) {
			Annotation annotation = new Annotation(new QName(addGeneNS,tag,"ag"), annotationValue);
			annotations.add(annotation);
		}
	}
	
	private static void addGeneStringAnnotationArray(List<Annotation> annotations, JSONObject object, String tag) throws SBOLValidationException {
		JSONArray annotationValues = (JSONArray)object.get(tag);
		for (Object annotationValue : annotationValues) {
			Annotation annotation = new Annotation(new QName(addGeneNS,tag,"ag"), (String)annotationValue);
			annotations.add(annotation);
		}
	}
	
	private static void createAddGeneStringAnnotation(TopLevel topLevel,JSONObject plasmid, String tag) throws SBOLValidationException {
		String annotation = plasmid.get(tag)!=null?plasmid.get(tag).toString():null;
		if (annotation!=null && !annotation.equals("")) {
			topLevel.createAnnotation(new QName(addGeneNS,tag,"ag"), annotation);
		}
	}
	
	private static void createAddGeneStringAnnotationArray(TopLevel topLevel,JSONObject object, String tag) throws SBOLValidationException {
		JSONArray annotationValues = (JSONArray)object.get(tag);
		for (Object annotationValue : annotationValues) {
			topLevel.createAnnotation(new QName(addGeneNS,tag,"ag"), (String)annotationValue);
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
	
	private static void createArticleAnnotations(ComponentDefinition cd, JSONObject plasmid) throws SBOLValidationException {
		JSONObject article = (JSONObject)plasmid.get("article");
		if (article != null) {
			if (article.get("pubmed_id")!=null) {
				String pubmedId = article.get("pubmed_id").toString();
				cd.createAnnotation(new QName(oboNS,"OBI_0001617","obo"), pubmedId);
			}
			if (article.get("doi")!=null) {
				String doi = article.get("doi").toString();
				cd.createAnnotation(new QName(oboNS,"OBI_0002110","obo"), doi);
			}
			if (article.get("id")!=null) {
				String id = article.get("id").toString();
				cd.createAnnotation(new QName(addGeneNS,"articleId","ag"), id);
			}
			if (article.get("url")!=null) {
				URI articleUrl = URI.create(article.get("url").toString());
				cd.createAnnotation(new QName(addGeneNS,"articleUrl","ag"), articleUrl);
			}
		}
	}
	
	private static void createCreatorAnnotations(ComponentDefinition cd, JSONObject plasmid) throws SBOLValidationException {
		JSONArray pis = (JSONArray)plasmid.get("pi");
		for (Object pi : pis) {
			cd.createAnnotation(new QName(dcNS,"creator","dc"), (String)pi);
		}
	}
	
	private static void createCloningAnnotations(ComponentDefinition cd, JSONObject plasmid) throws SBOLValidationException {
		JSONObject cloning = (JSONObject)plasmid.get("cloning");
		List<Annotation> annotations = new ArrayList<Annotation>();
		addGeneStringAnnotation(annotations,cloning,"backbone");
		addGeneStringAnnotation(annotations,cloning,"backbone_mutation");
		addGeneStringAnnotation(annotations,cloning,"backbone_origin");
		addGeneStringAnnotation(annotations,cloning,"backbone_size");
		addGeneStringAnnotation(annotations,cloning,"promoter");
		addGeneStringAnnotation(annotations,cloning,"sequencing_primer_3");
		addGeneStringAnnotation(annotations,cloning,"sequencing_primer_5");
		addGeneStringAnnotationArray(annotations,cloning,"vector_types");
		cd.createAnnotation(new QName(addGeneNS,"cloning","ag"), new QName(addGeneNS,"Cloning","ag"), 
				URI.create(cd.getPersistentIdentity().toString()+"/cloning"), annotations);
	}
	
	private static void createTagAnnotation(List<Annotation> tagAnnotations, JSONObject object, 
			String nestedUriPrefix) throws SBOLValidationException {
		JSONArray tags = (JSONArray)object.get("tags");
		int i = 0;
		for (Object tag : tags) {
			List<Annotation> annotations = new ArrayList<Annotation>();
			addGeneStringAnnotation(annotations,(JSONObject)tag,"location");
			addGeneStringAnnotation(annotations,(JSONObject)tag,"tag");
			Annotation annotation = new Annotation(new QName(addGeneNS,"tagElement","ag"), new QName(addGeneNS,"TagElement","ag"), 
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
		addGeneStringAnnotation(annotations,(JSONObject)cloning,"clone_method");
		addGeneStringAnnotation(annotations,(JSONObject)cloning,"cloning_site_3");
		addGeneStringAnnotation(annotations,(JSONObject)cloning,"cloning_site_5");
		addGeneStringAnnotation(annotations,(JSONObject)cloning,"promoter");
		addGeneStringAnnotation(annotations,(JSONObject)cloning,"sequencing_primer_3");
		addGeneStringAnnotation(annotations,(JSONObject)cloning,"sequencing_primer_5");
		addGeneStringAnnotation(annotations,(JSONObject)cloning,"site_3_destroyed");
		addGeneStringAnnotation(annotations,(JSONObject)cloning,"site_5_destroyed");
		Annotation annotation = new Annotation(new QName(addGeneNS,"cloning","ag"), new QName(addGeneNS,"Cloning","ag"), 
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
			addGeneStringAnnotation(annotations,(JSONObject)entrezGene,"aliases");
			addGeneStringAnnotation(annotations,(JSONObject)entrezGene,"gene");
			addGeneStringAnnotation(annotations,(JSONObject)entrezGene,"id");
			Annotation annotation = new Annotation(new QName(addGeneNS,"entrezGene","ag"), new QName(addGeneNS,"EntrezGene","ag"), 
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
					Annotation annotation = new Annotation(new QName(addGeneNS,"species","ag"), species.toString());
					insertAnnotations.add(annotation);
				}
			}
		}
	}
	
	private static void createInsertAnnotations(ComponentDefinition cd, JSONObject plasmid) throws SBOLValidationException {
		JSONArray inserts = (JSONArray)plasmid.get("inserts");
		int i = 0;
		for (Object insert : inserts) {
			List<Annotation> annotations = new ArrayList<Annotation>();
			addGeneStringAnnotationArray(annotations,(JSONObject)insert,"alt_names");
			createInsertCloningAnnotations(annotations,(JSONObject)insert,i,cd.getPersistentIdentity().toString());
			createInsertEntrezGeneAnnotations(annotations,(JSONObject)insert,i,cd.getPersistentIdentity().toString());
			createInsertSpeciesAnnotations(annotations,(JSONObject)insert);
			addGeneStringAnnotationArray(annotations,(JSONObject)insert,"genbank_ids");
			addGeneStringAnnotation(annotations,(JSONObject)insert,"mutation");
			addGeneStringAnnotation(annotations,(JSONObject)insert,"name");
			addGeneStringAnnotation(annotations,(JSONObject)insert,"shRNA_sequence");
			addGeneStringAnnotation(annotations,(JSONObject)insert,"size");
			createTagAnnotation(annotations,(JSONObject)insert,cd.getPersistentIdentity().toString()+"/insert"+i);
			cd.createAnnotation(new QName(addGeneNS,"insert","ag"), new QName(addGeneNS,"Insert","ag"), 
					URI.create(cd.getPersistentIdentity().toString()+"/insert"+i), annotations);
			i++;
		}
	}
	
	public static void main( String[] args ) throws FileNotFoundException, IOException, ParseException 
    {
		// Create an SBOLDocument
		SBOLDocument document = new SBOLDocument(); 
		document.setDefaultURIprefix(uriPrefix); 
		document.setComplete(true); 
		document.setCreateDefaults(true);
		SynBioHubFrontend sbh = new SynBioHubFrontend("https://synbiohub.utah.edu");
		URI activity = null;
//		SynBioHubFrontend sbh = new SynBioHubFrontend("http://localhost:7777","https://synbiohub.org");
		int start = 14372;
		
		try {
			// Create an Activity
			GenericTopLevel genericTopLevel = document.createGenericTopLevel("addgene2sbol", version, 
					new QName(provNS, "Activity", "prov"));
			activity = genericTopLevel.getIdentity();
			genericTopLevel.setName("AddGene to SBOL conversion");
			genericTopLevel.setDescription("Conversion of the Addgene plasmids and metadata to SBOL2");
			genericTopLevel.createAnnotation(new QName(dcNS,"creator","dc"), "Chris J. Myers");
			TimeZone tz = TimeZone.getTimeZone("UTC");
			DateFormat df = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss.SSS'Z'");
			df.setTimeZone(tz);
			createdDate = df.format(new Date());
			genericTopLevel.createAnnotation(new QName(provNS,"endedAtTime","prov"), createdDate);
			activityURI = genericTopLevel.getIdentity();

			// Create collection
			System.out.println(args[0]);
			System.out.println(args[1]);
			System.out.println(args[2]);
			sbh.login(args[0], args[1]);
			if (start==0)
				sbh.createCollection("AddGenePlasmids", "1", "AddGene Plasmids", "These are the AddGene plasmids", "", true, document);
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
		for (Object p : plasmids) {
			i++;
			if (i < start) continue;
			document = new SBOLDocument(); 
			document.setDefaultURIprefix(uriPrefix); 
			//document.setComplete(true); 
			document.setCreateDefaults(true);
			JSONObject plasmid = (JSONObject) p;
			String displayId = "addGene"+plasmid.get("id");
			try {
				ComponentDefinition cd = document.createComponentDefinition(displayId, version, ComponentDefinition.DNA);
				cd.addWasGeneratedBy(URI.create("https://synbiohub.utah.edu/user/myers/AddGenePlasmids/addgene2sbol/1"));
				cd.addWasDerivedFrom(URI.create(plasmid.get("url").toString()));
				cd.setName(plasmid.get("name").toString());
				cd.addRole(URI.create(so + "SO:0000155"));
				createCreatorAnnotations(cd,plasmid);
				createArticleAnnotations(cd,plasmid);
				createCloningAnnotations(cd,plasmid);
				createInsertAnnotations(cd,plasmid);
				List<Annotation> tagAnnotations = new ArrayList<Annotation>();
				createTagAnnotation(tagAnnotations,plasmid,cd.getPersistentIdentity().toString());
				if (tagAnnotations.size() > 0) {
					cd.createAnnotation(new QName(addGeneNS,"tags","ag"), new QName(addGeneNS,"Tags","ag"), 
							URI.create(cd.getPersistentIdentity()+"/tag"), tagAnnotations);
				}
				createAddGeneStringAnnotation(cd,plasmid,"bacterial_resistance");
				createAddGeneStringAnnotation(cd,plasmid,"growth_notes");
				createAddGeneStringAnnotation(cd,plasmid,"growth_strain");
				createAddGeneStringAnnotation(cd,plasmid,"growth_temp");
				createAddGeneStringAnnotation(cd,plasmid,"origin");
				createAddGeneStringAnnotation(cd,plasmid,"plasmid_copy");
				createAddGeneStringAnnotationArray(cd,plasmid,"resistance_markers");
				createAddGeneStringAnnotation(cd,plasmid,"terms");
				JSONObject sequences = (JSONObject) plasmid.get("sequences");
				createSequence(document,cd,sequences,"public_addgene_full_sequences","_addgene_full_seq");
				createSequence(document,cd,sequences,"public_addgene_partial_sequences","_addgene_partial_seq");
				createSequence(document,cd,sequences,"public_user_full_sequences","_user_full_seq");
				createSequence(document,cd,sequences,"public_user_partial_sequences","_user_partial_seq");
				SBOLValidate.validateSBOL(document,false,true,false);
				if (SBOLValidate.getNumErrors()>0) {
					failure++;
					System.out.println(i + " out of " + size + ":"+displayId+" FAILURE "+ failure);
					System.err.println(i + " out of " + size + ":"+displayId+" FAILURE "+ failure);
					for (String error : SBOLValidate.getErrors()) {
						System.err.println(error);
					}
				} else {   
					// Upload to SynBioHub
					sbh.addToCollection(URI.create("https://synbiohub.utah.edu/user/myers/AddGenePlasmids/AddGenePlasmids_collection/1"), true, document);
					success++;
		        	System.out.println(i + " out of " + size + ":"+displayId+" SUCCESS "+ success);
				}
			} catch (Exception e) {
				failure++;
	        	System.out.println(i + " out of " + size + ":"+displayId+" FAILURE "+ failure);
	        	System.err.println(i + " out of " + size + ":"+displayId+" FAILURE "+ failure);
	        	e.printStackTrace(System.err);
			}
		}
		//document.write(System.out);
        
        // Validate
    }
}
