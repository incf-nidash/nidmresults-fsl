'''Python implemetation of the export of FSL results into NI-DM

@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013
'''

from HTMLParser import HTMLParser
from htmlentitydefs import name2codepoint
import re
from prov.model import ProvBundle, ProvRecord, ProvExceptionCannotUnifyAttribute, graph, ProvEntity
import prov.model.graph

class FSL_NIDM():

    def __init__(self, *args, **kwargs):
        self.reportFile = []
        self.zstatFile = []
        self.create_basis_fsl_prov()

    def create_basis_fsl_prov(self):
        g = ProvBundle()
        g.add_namespace("neurolex", "http://neurolex.org/wiki/")
        g.add_namespace("spm", "http://www.fil.ion.ucl.ac.uk/spm/ns/")
        g.add_namespace("nidm", "http://nidm.nidash.org/")
        g.add_namespace("niiri", "http://iri.nidash.org/")
        g.add_namespace("crypto", "http://www.w3.org/2000/10/swap/crypto#")

        g.entity('niiri:statistical_map_id' ,
            other_attributes=(  ('prov:type', 'nidm:statisticalMap'), 
                                ('prov:label', "Statistical Map: TODO") ,
                                ('prov:contrastName', "TODO") ,
                                ('prov:location', "file:///path/to/TODO.img"),
                                ('prov:fileName', "TODO.img"),
                                ('prov:statisticType', 'nidm:tStatisticTODO'),
                                ('prov:errorDegreesOfFreedom', 'TODO'),
                                ('prov:effectDegreesOfFreedom', 'TODO')
                                ) )
        g.entity('niiri:contrast_map_id', other_attributes=(  ('prov:type', 'nidm:contrastMap'), ('prov:type', 'nidm:contrastMap2')))
        g.wasDerivedFrom('niiri:statistical_map_id', 'niiri:contrast_map_id')

        g.entity('niiri:contrast_standard_error_map_id', other_attributes=( ('prov:type' , 'nidm:contrastStandardErrorMap'), ('prov:type' , 'nidm:contrastStandardErrorMap2')))

        g.wasDerivedFrom('niiri:statistical_map_id', 'niiri:contrast_standard_error_map_id')
        g.entity('niiri:residual_mean_squares_map_id', other_attributes=( ('prov:type','nidm:residualMeanSquaresMap',), ('prov:location',"file:///path/to/ResMS.img" )))
        g.wasDerivedFrom('niiri:contrast_standard_error_map_id', 'niiri:residual_mean_squares_map_id')
        g.entity('niiri:design_matrix_id', other_attributes=( ('prov:type','nidm:designMatrix',), ('prov:location', "file:///path/to/design_matrix.csv")))
        g.entity('niiri:contrast_id', other_attributes=( ('prov:type', 'nidm:contrast'), ('nidm:contrastName',"listening &gt; rest")))
        g.wasDerivedFrom('niiri:contrast_map_id', 'niiri:design_matrix_id')
        g.wasDerivedFrom('niiri:contrast_map_id', 'niiri:contrast_id')
        g.wasDerivedFrom('niiri:contrast_standard_error_map_id', 'niiri:design_matrix_id')
        g.wasDerivedFrom('niiri:contrast_standard_error_map_id', 'niiri:contrast_id')
        g.wasDerivedFrom('niiri:residual_mean_squares_map_id', 'niiri:design_matrix_id')
        g.entity('niiri:search_space_id',other_attributes=( ('prov:type', 'nidm:mask',), ('prov:location', "file:///path/to/mask.img")))
        
        g.activity('niiri:inference_id', other_attributes=( ('prov:type', 'spm:inference'), ('prov:label' , "Inference")))
        g.used('niiri:inference_id', 'niiri:statistical_map_id')
        g.used('niiri:inference_id', 'niiri:search_space_id')
        g.used('niiri:inference_id', 'niiri:height_threshold_id')
        g.used('niiri:inference_id', 'niiri:extent_threshold_id')
        
        
        g.entity('niiri:stat_image_properties_id', other_attributes=( ('prov:type', 'spm:statisticImageProperties'), ('spm:expectedNumberOfVoxelsPerCluster', "0.553331387916112")))
        g.wasGeneratedBy('niiri:stat_image_properties_id', 'niiri:inference_id')
        g.entity('niiri:excursion_set_id', other_attributes=( 
            ('prov:type', 'fsl:excursionSet'), 
            ('prov:location',"file:///path/to/thresh_zstat1.nii.gz"),
            ('nidm:voxelToWorldMapping', "[[1,2,3,4], [1,2,3,4], [1,2,3,4], [0,0,0,1]]"),
            ('nidm:numberOfDimensions', "3"),
            ('nidm:dimensions', "[53 63 46]"),
            ('nidm:voxelUnits', "['mm' 'mm' 'mm']"),
            ('nidm:voxelSize', "[3 3 3]"),
            ))
        g.wasGeneratedBy('niiri:excursion_set_id', 'niiri:inference_id')
        self.provBundle = g

    def add_report_file(self, myReportFile):
        self.reportFile = myReportFile
        parser = MyFSLReportParser();
        file = open(myReportFile, 'r')
        parser.feed(file.read());

        heightThreshAllFields = {
            'prov:type': 'nidm:heightThreshold', 'prov:label': "Height Threshold: TODO",
            'prov:userSpecifiedThresholdType': "TODO" , 'prov:value': parser.get_voxelThreshValue(),
            'nidm:pValueUncorrected': parser.get_voxelPUncorr(), 'fsl:pValueGRF': parser.get_voxelPCorr()
            }
        self.provBundle.entity('niiri:height_threshold_id', other_attributes=dict((k,v) for k,v in heightThreshAllFields.iteritems() if v is not None))
        exentThreshAllFields = {
            'prov:type': 'nidm:extentThreshold', 'prov:label': "Extent Threshold: TODO", 'nidm:clusterSizeInVoxels': parser.get_extentValue(),
            'nidm:pValueUncorrected': parser.get_extentPUncorr(), 'fsl:pValueGRF': parser.get_extentPCorr()
        }
        self.provBundle.entity('niiri:extent_threshold_id', other_attributes=dict((k,v) for k,v in exentThreshAllFields.iteritems() if v is not None))

        self.provBundle.agent('niiri:software_id', other_attributes=( 
            ('prov:type', 'nidm:fsl'), 
            ('prov:type','prov:SoftwareAgent'),
            ('prov:label','FSL'),
            ('nidm:softwareVersion','FSL') ))

        self.provBundle.wasAssociatedWith('niiri:inference_id', 'niiri:software_id')


    def add_zstat_file(self, myZstatFile, myStdZstatFile):
        self.zstatFile = myZstatFile
        clusterParser = MyFSLClusterParser();
        fileVox = open(myZstatFile, 'r')
        clusterParser.feed(fileVox.read());
        clusterParserInStdSpace = MyFSLClusterParser();
        fileUnits = open(myStdZstatFile, 'r')
        clusterParserInStdSpace.feed(fileUnits.read());
        
        clusters = clusterParser.get_clusters()
        clustersStd = clusterParserInStdSpace.get_clusters()
        clusIdx = -1
        if clusters is not None:
            for cluster in clusters:               
                clusIdx = clusIdx + 1
                self.provBundle.entity('niiri:cluster_'+str(cluster.get_id()), other_attributes=( 
                             ('prov:type' , 'spm:clusterLevelStatistic'), 
                             ('prov:label', "Cluster Level Statistic: "+str(cluster.get_id())),
                             ('fsl:clusterSizeInVoxels', str(cluster.get_SizeInVoxels())),
                             ('fsl:pGRF', str(cluster.get_pGRF()) )))
                self.provBundle.wasDerivedFrom('niiri:cluster_'+str(cluster.get_id()), 'niiri:excursion_set_id')

                self.provBundle.entity('niiri:coordinate_'+str(cluster.get_id())+"00", other_attributes=( 
                    ('prov:type' , 'prov:location'), 
                    ('prov:type' , 'nidm:coordinate'),
                    ('nidm:coordinateSystem' , 'nidm:mni'),
                    ('nidm:coordinate1' , str(cluster.get_COG1())),
                    ('nidm:coordinate2' , str(cluster.get_COG2())),
                    ('nidm:coordinate3' , str(cluster.get_COG3())),
                    ('nidm:coordinateInUnits1' , str(clustersStd[clusIdx].get_COG1())),
                    ('nidm:coordinateInUnits2' , str(clustersStd[clusIdx].get_COG2())),
                    ('nidm:coordinateInUnits3' , str(clustersStd[clusIdx].get_COG3()))
                    ))
                self.provBundle.entity('niiri:centerOfGravity_'+str(cluster.get_id()), other_attributes=( 
                             ('prov:type' , 'fsl:centerOfGravity'), 
                             ('prov:label', "Center of Gravity: "+str(cluster.get_id())),
                             ('prov:atLocation' , 'niiri:coordinate_'+str(cluster.get_id())+"00"))   )
                self.provBundle.wasDerivedFrom('niiri:centerOfGravity_'+str(cluster.get_id()), 'niiri:cluster_'+str(cluster.get_id()))
        peaks = clusterParser.get_peaks()
        peaksStd = clusterParserInStdSpace.get_peaks()
        if peaks is not None:
            peakIndex = 1;
            for peak in peaks:      
                self.provBundle.entity('niiri:coordinate_'+str(peakIndex), other_attributes=( 
                    ('prov:type' , 'prov:location'), 
                    ('prov:type' , 'nidm:coordinate'),
                    ('prov:label' , "Coordinate "+str(peakIndex)),
                    ('nidm:coordinate1' , peak.get_x()),
                    ('nidm:coordinate2' , peak.get_y()),
                    ('nidm:coordinate3' , peak.get_z()),
                    ('nidm:coordinateInUnits1' , str(peaksStd[peakIndex-1].get_x())),
                    ('nidm:coordinateInUnits2' , str(peaksStd[peakIndex-1].get_y())),
                    ('nidm:coordinateInUnits3' , str(peaksStd[peakIndex-1].get_z()))
                    ))
                self.provBundle.entity('niiri:peak_'+str(peakIndex), other_attributes=( 
                    ('prov:type' , 'spm:peakLevelStatistic'), 
                    ('nidm:equivalentZStatistic', str(peak.get_equivZStat())),
                    ('prov:atLocation' , 'niiri:coordinate_'+str(peakIndex)))         )
                self.provBundle.wasDerivedFrom('niiri:peak_'+str(peakIndex), 'niiri:cluster_'+str(peak.get_cluster_id()))
                peakIndex = peakIndex + 1
        

    def save_prov_to_files(self):
        jsondata = self.provBundle.get_provjson(indent=4)
        JSONfile = open('./FSL_example.json', 'w');
        JSONfile.write(jsondata)
        PROVNfile = open('./FSL_example.provn', 'w');
        PROVNfile.write(self.provBundle.get_provn(4))

        dot = graph.prov_to_dot(self.provBundle)
        # , use_labels=True)
        dot.set_dpi(120)
        dot.write_png('./FSL_example.png')

'''HTML parser for report files

'''
class MyFSLReportParser(HTMLParser):

    def __init__(self, *args, **kwargs):
        HTMLParser.__init__(self, *args, **kwargs)
        self.descriptions = []
        self.inside_a_element = 0
        self.hyperlinks = []
        self.foundIntro = False;
        self.featVersion = ''
        self.pValue = []
        self.threshType = ''

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for name, value in attrs:
                if name == "href":
                    self.hyperlinks.append(value)
                    self.inside_a_element = 1

    def handle_endtag(self, tag):
        if tag == 'a':
            self.inside_a_element = 0
    def handle_data(self, data):
        if self.inside_a_element:
            self.descriptions.append(data)
        elif not self.foundIntro:
            # Look for p-value, type of thresholding and feat version in introductory text
            patternVoxelThresh = re.compile(r'.*Version (?P<featversion>\d+\.\d+),.* thresholded using (?P<threshtype>.*) thresholding .* P=(?P<pvalue>\d+\.\d+)')

            extractedData = patternVoxelThresh.search(data) 
            
            if extractedData is not None:
                self.featVersion = extractedData.group('featversion')
                self.voxelThreshValue = None;
                self.voxelPCorr = extractedData.group('pvalue')
                self.voxelPUncorr = None
                self.extentValue = 0;
                self.extentPCorr = 1
                self.extentPUncorr = 1
                # self.threshType = extractedData.group('threshtype')
                self.foundIntro = True;
            else:
                patternClusterThresh = re.compile(r'.*Version (?P<featversion>\d+\.\d+),.* thresholded using (?P<threshtype>.*) determined by Z\>(?P<zvalue>\d+\.\d+) and a .* P=(?P<pvalue>\d+\.\d+) .*')
                extractedData = patternClusterThresh.search(data) 

                if extractedData is not None:
                    self.featVersion = extractedData.group('featversion')
                    self.voxelThreshValue = extractedData.group('zvalue')
                    self.voxelPCorr = None
                    self.voxelPUncorr = None
                    self.extentValue = None;
                    self.extentPCorr = extractedData.group('pvalue');
                    self.extentPUncorr = None
                    # self.threshType = extractedData.group('threshtype')
                    self.foundIntro = True;

    def get_threshold_p_value(self):
        return self.pValue

    def get_voxelThreshValue(self):
        return self.voxelThreshValue

    def get_voxelPCorr(self):
        return self.voxelPCorr

    def get_voxelPUncorr(self):
        return self.voxelPUncorr

    def get_extentValue(self):
        return self.extentValue

    def get_extentPCorr(self):
        return self.extentPCorr

    def get_extentPUncorr(self):
        return self.extentPUncorr

'''HTML parser for cluster files

'''
class MyFSLClusterParser(HTMLParser):
    def __init__(self, *args, **kwargs):
        HTMLParser.__init__(self, *args, **kwargs)
        self.descriptions = None
        self.inside_a_element = 0
        self.hyperlinks = None
        self.inside_h1 = False;
        self.inside_cluster_list = False;
        self.inside_peak_list = False;
        self.in_cluster_table = False;
        self.in_peak_table = False;
        self.row_number = 0;
        self.clusters = []
        self.clusterAppended = False
        self.peaks = []
        self.peakAppended = False
        # self.peakIndex = 1

    def handle_starttag(self, tag, attrs):
        if tag == "h1":
            self.inside_h1 = True;
            self.inside_cluster_list = False;            
        elif tag == "table":
            if self.inside_cluster_list:         
                self.in_cluster_table = True;
                self.row_number = 0;
            elif self.inside_peak_list:         
                self.in_peak_table = True;
                self.row_number = 0;
        elif tag == "tr":
            self.row_number = self.row_number + 1;
            self.col_number = 0;
        elif tag == "td":
            self.col_number = self.col_number + 1;

    def handle_endtag(self, tag):
        if tag == "h1":    
            self.inside_h1 = False;
        elif tag == "table":
            self.in_cluster_table = False;
        

    def handle_data(self, data):
        if self.inside_h1:
            if data == "Cluster List":
                self.inside_cluster_list = True;
            elif data == "Local Maxima":
                self.inside_peak_list = True;

        if self.in_cluster_table:
            if self.row_number > 1:
                if self.col_number == 1:
                    self.cluster = Cluster(data);
                    self.clusterAppended = False
                elif self.col_number == 2:
                    self.cluster.sizeInVoxels(data)
                elif self.col_number == 3:
                    self.cluster.set_pGRF(data)
                elif self.col_number == 9:
                    self.cluster.set_COG1(data)
                elif self.col_number == 10:
                    self.cluster.set_COG2(data)
                elif self.col_number == 11:
                    self.cluster.set_COG3(data)
                elif self.col_number == 16:
                    if not self.clusterAppended:
                        self.clusters.append(self.cluster)
                        self.clusterAppended = True
        

        if self.in_peak_table:
            if self.row_number > 1:
                if self.col_number == 1:
                    self.peak = Peak(data);
                    self.peakAppended = False
                elif self.col_number == 2:
                    self.peak.set_equivZStat(data)
                elif self.col_number == 3:
                    self.peak.set_x(data)
                elif self.col_number == 4:
                    self.peak.set_y(data)
                elif self.col_number == 5:
                    if not self.peakAppended:
                        self.peak.set_z(data)
                        self.peaks.append(self.peak)
                        self.peakAppended = True

                        
    def get_clusters(self):
        return self.clusters

    def get_peaks(self):
        return self.peaks        

class Peak():
    def __init__(self, clusterId, *args, **kwargs):
        self.cluster_id = clusterId
        self.equivZStat = None
        self.x = None
        self.y = None
        self.z = None

    def set_equivZStat(self,value):
        self.equivZStat = value

    def get_cluster_id(self):
        return self.cluster_id    

    def get_equivZStat(self):
        return self.equivZStat    

    def set_x(self,value):
        self.x = value
    def set_y(self,value):
        self.y = value
    def set_z(self,value):
        self.z = value

    def get_x(self):
        return self.x  
    def get_y(self):
        return self.y  
    def get_z(self):
        return self.z                  

class Cluster():
    def __init__(self, clusterId, *args, **kwargs):
        self.cluster_id = clusterId
        self.clusterSizeInVoxels = None
        self.equivZStat = None
        self.pGRF = None
        self.COG1 = None
        self.COG2 = None
        self.COG3 = None

    def sizeInVoxels(self,value):
        self.clusterSizeInVoxels = value

    

    def set_COG1(self,value):
        self.COG1 = value   
    def set_COG2(self,value):
        self.COG2 = value   
    def set_COG3(self,value):
        self.COG3 = value   

    def get_COG1(self):
        return self.COG1
    def get_COG2(self):
        return self.COG2 
    def get_COG3(self):
        return self.COG3            

    def set_pGRF(self,value):
        self.pGRF = value        

    def get_id(self):
        return self.cluster_id

    def get_SizeInVoxels(self):
        return self.clusterSizeInVoxels



    def get_pGRF(self):
        return self.pGRF



