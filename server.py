#!/usr/bin/env python

import os, uuid, urllib, re, sys, datetime, webbrowser
import tornado.httpserver
import tornado.ioloop
import tornado.options
import tornado.web
import tornado.websocket
from pymongo import MongoClient
from bson.objectid import ObjectId

mongo_client = None

class NCBI:

    def __init__(self):
        self._eutils_base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    def get_entry(self, entry_id):
        response = urllib.urlopen("%sefetch.fcgi?db=nucleotide&id=%s&rettype=gbwithparts"%(self._eutils_base_url, entry_id))
        content = str(response.read())
        response.close()
        return content

    def parse_entry(self, entry_content):
        pieces_of_seq=[]
        start_of_sequence = False
        accession = None
        feature_type = None
        qualifer_type = None
        qualifier_content = None
        qualifiers = []
        genomic_strand = '+'
        genomic_positions = None
        features = []
        organism = None
        inOrganism = False
        lineage = ""
        location = None
        lines = entry_content.strip().split('\n')
        if not lines[-1].strip() == '//':
            raise Exception("Uncomplete file")
        for line in lines:
            tokens = re.split('\s+', line)
            if line.startswith('ACCESSION'):
                accession = re.split('\s+', line)[1]
            elif line.strip().startswith('ORGANISM'):
                organism = line.strip().split('ORGANISM')[1].strip()
                inOrganism = True
            elif line.strip().startswith('REFERENCE'):
                inOrganism = False
            elif line.startswith('ORIGIN'): #the genomic sequence
                start_of_sequence = True
                #we store the last feature
                #the last
                if feature_type and not feature_type == "source" and not feature_type == "intron":
                    if location.startswith("complement(join("):#a joined location on the Crick strand
                        genomic_strand = '-'
                        ends = location.split('complement(join(')[1][:-2].replace("join(", "").replace(')','').replace(',','..').split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]
                        if feature_type == 'CDS':
                            for i in range(1, len(ends)-2, 2):
                                intron = {
                                    'type': 'intron',
                                    'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                    'genomicStrand': genomic_strand,
                                }
                                features.append(intron)
                    elif location.startswith("join(complement("):#a joined location on the Crick strand
                        genomic_strand = '-'
                        ends = location.split('join(complement(')[1][:-2].replace("complement(", "").replace(')', '').replace(',','..').split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]
                        if feature_type == 'CDS':
                            for i in range(1, len(ends)-2, 2):
                                intron = {
                                    'type': 'intron',
                                    'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                    'genomicStrand': genomic_strand,
                                }
                                features.append(intron)
                    elif location.startswith("complement(order("):
                        genomic_strand = '-'
                        ends = location.split('complement(order(')[1][:-2].replace(',','..').split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]
                    elif location.startswith("order("):
                        ends = location.split('order(')[1][:-1].replace(',','..').split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]
                    elif location.startswith("complement("): #a location on the Crick strand
                        genomic_strand = '-'
                        ends = location.split('complement(')[1][:-1].split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]
                    elif location.startswith("join("): #a joined location
                        ends = location.split('join(')[1][:-1].replace(',','..').split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]
                        if feature_type == 'CDS':
                            for i in range(1, len(ends)-2, 2):
                                intron = {
                                    'type': 'intron',
                                    'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                    'genomicStrand': genomic_strand,
                                }
                                features.append(intron)
                    else: #a regular location
                        ends = location.split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]

                    feature = {
                        'type': feature_type,
                        'genomicPositions': genomic_positions,
                        'genomicStrand': genomic_strand,
                    }
                    if qualifer_type and qualifier_content:
                        if qualifer_type == 'translation':
                            qualifier_content = qualifier_content.replace(" ","")
                        qualifiers.append({
                            "type": qualifer_type,
                            "content": qualifier_content
                             })
                    for qualifier in qualifiers:
                        feature[qualifier['type']] = qualifier['content']
                    features.append(feature)
            elif len(tokens) == 3 and re.findall('\.\.>?[0-9]+', tokens[2]):
                #new feature
                #we store the previous one (if any)
                if feature_type and not feature_type == "source" and not feature_type == "intron":
                    if location.startswith("complement(join("):#a joined location on the Crick strand
                        genomic_strand = '-'
                        ends = location.split('complement(join(')[1][:-2].replace("join(", "").replace(')','').replace(',','..').split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]
                        if feature_type == 'CDS':
                            for i in range(1, len(ends)-2, 2):
                                intron = {
                                    'type': 'intron',
                                    'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                    'genomicStrand': genomic_strand,
                                }
                                features.append(intron)
                    elif location.startswith("join(complement("):#a joined location on the Crick strand
                        genomic_strand = '-'
                        ends = location.split('join(complement(')[1][:-2].replace("complement(", "").replace(')', '').replace(',','..').split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]
                        if feature_type == 'CDS':
                            for i in range(1, len(ends)-2, 2):
                                intron = {
                                    'type': 'intron',
                                    'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                    'genomicStrand': genomic_strand,
                                }
                                features.append(intron)
                    elif location.startswith("complement(order("):
                        genomic_strand = '-'
                        ends = location.split('complement(order(')[1][:-2].replace(',','..').split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]
                    elif location.startswith("order("):
                        ends = location.split('order(')[1][:-1].replace(',','..').split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]
                    elif location.startswith("complement("): #a location on the Crick strand
                        genomic_strand = '-'
                        ends = location.split('complement(')[1][:-1].split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]
                    elif location.startswith("join("): #a joined location
                        ends = location.split('join(')[1][:-1].replace(',','..').split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]
                        if feature_type == 'CDS':
                            for i in range(1, len(ends)-2, 2):
                                intron = {
                                    'type': 'intron',
                                    'genomicPositions': [ends[i]+1, ends[i+1]-1],
                                    'genomicStrand': genomic_strand,
                                }
                                features.append(intron)
                    else: #a regular location
                        ends = location.split('..')
                        ends = map(lambda end: int(end.replace('>','').replace('<','')), ends)
                        genomic_positions = [min(ends), max(ends)]

                    feature = {
                        'type': feature_type,
                        'genomicPositions': genomic_positions,
                        'genomicStrand': genomic_strand,
                    }

                    if qualifer_type and qualifier_content:
                        if qualifer_type == 'translation':
                            qualifier_content = qualifier_content.replace(" ","")
                        qualifiers.append({
                            "type": qualifer_type,
                            "content": qualifier_content
                             })
                    for qualifier in qualifiers:
                        feature[qualifier['type']] = qualifier['content']
                    features.append(feature)
                feature_type = None
                genomic_strand = '+'
                genomic_positions = None
                qualifer_type = None
                qualifier_content = None
                qualifiers = []
                feature_type = tokens[1].strip()
                location = tokens[2].strip()
            elif not qualifer_type and not qualifier_content and len(tokens) == 2 and re.findall('\.\.',tokens[1]): #still the content of the current location
                location += tokens[1].strip()
            elif re.findall('^\s+/.+=', line): # a new qualifier /bla_bla=
                if qualifer_type and qualifier_content:
                    if qualifer_type == 'translation':
                        qualifier_content = qualifier_content.replace(" ","")
                    qualifiers.append({
                        "type": qualifer_type,
                        "content": qualifier_content
                         })
                qualifer_type = line.strip()[1:].split('=')[0].strip()[0:]
                qualifier_content = line.strip()[1:].split('=')[1].strip().replace('"','')
            elif re.findall('^\s+/.+', line): # a qualifier like /manual => ignore
                pass
            elif not start_of_sequence and qualifer_type and qualifier_content : #still the content of the current qualifier
                qualifier_content += " "+line.strip().replace('"','')
            elif line.startswith('//'): #end of the genomic sequence
                start_of_sequence = False
            elif start_of_sequence:
                pieces_of_seq.append(''.join(re.split('\s+',line.strip())[1:]).upper())
            elif inOrganism:
                lineage += " "+line.strip()

        return organism, ''.join(pieces_of_seq), features

class IndexHandler(tornado.web.RequestHandler):
    def get(self):
        self.render('index.html')

class WebSocketHandler(tornado.websocket.WebSocketHandler):
    clients = {}

    def open(self, *args):
        self.id = uuid.uuid4()
        self.clients[self.id] = {'id':self.id}
        print "New client connected"

    def on_message(self, message):
        import json
        message = json.loads(message)
        if message['header'] == 'get available projects':
            databases_names  = mongo_client.database_names()
            answer = {'header': 'got available projects'}
            answer['projects'] = databases_names
            self.write_message(answer, binary = False)
        
        elif message['header'] == 'create project':
            project_name = message['name'].replace(' ','_')
            ncbi = NCBI()
            organism, sequence, features = ncbi.parse_entry(ncbi.get_entry(message['ncbi_id']))

            self.clients[self.id]['db'] = mongo_client[project_name]

            genome_description = {
                '_id': str(ObjectId()),
                'name': message['ncbi_id'],
                'sequence': sequence,
                'source': 'db:ncbi:%s'%message['ncbi_id'],
                'organism': organism
            }
        
            self.clients[self.id]['db']['genomes'].insert(genome_description)

            annotations = []

            for feature in features:
                annotation = {
                    '_id': str(ObjectId()),
                    'class': feature['type'],
                    'source': 'db:ncbi:%s'%message['ncbi_id'],
                    'organism': organism,
                    'genomicPositions': feature['genomicPositions'],
                    'genomicStrand': feature['genomicStrand'],
                    'genome': "%s@genomes"%genome_description['_id'],
                    'genomeName': message['ncbi_id'],
                }

                for key in feature:
                    if not key in ['type', 'genomicPositions', 'genomicStrand'] and feature.get(key,None):
                        annotation[key] = feature[key]

                annotations.append(annotation)

                self.clients[self.id]['db']['annotations'].insert(annotation)

            answer = {
                'header': 'project created',
                'project': project_name,
                }
            answer['annotations'] = annotations
            self.write_message(answer, binary = False)

        elif message['header'] == 'add data':
            ncbi = NCBI()
            organism, sequence, features = ncbi.parse_entry(ncbi.get_entry(message['ncbi_id']))

            self.clients[self.id]['db'] = mongo_client[message['project']]

            genome_description = {
                '_id': str(ObjectId()),
                'name': message['ncbi_id'],
                'sequence': sequence,
                'source': 'db:ncbi:%s'%message['ncbi_id'],
                'organism': organism
            }
        
            self.clients[self.id]['db']['genomes'].insert(genome_description)

            annotations = []

            for feature in features:
                annotation = {
                    '_id': str(ObjectId()),
                    'class': feature['type'],
                    'source': 'db:ncbi:%s'%message['ncbi_id'],
                    'organism': organism,
                    'genomicPositions': feature['genomicPositions'],
                    'genomicStrand': feature['genomicStrand'],
                    'genome': "%s@genomes"%genome_description['_id'],
                    'genomeName': message['ncbi_id'],
                }

                for key in feature:
                    if not key in ['type', 'genomicPositions', 'genomicStrand'] and feature.get(key,None):
                        annotation[key] = feature[key]

                annotations.append(annotation)

                self.clients[self.id]['db']['annotations'].insert(annotation)

            answer = {
                'header': 'data added'
                }
            answer['annotations'] = annotations
            self.write_message(answer, binary = False)    
                
        elif message['header'] == 'remove project':
            mongo_client.drop_database(message['project'])
            answer = {
                'header': 'project removed'
                }
            self.write_message(answer, binary = False) 
               
        elif message['header'] == 'get all annotations':
            project = message['project']
            self.clients[self.id]['db'] = mongo_client[project]
            annotations = []
            for annotation in self.clients[self.id]['db']['annotations'].find():
                annotations.append(annotation)
            answer = {'header': 'got all annotations'}
            answer['annotations'] = annotations
            self.write_message(answer, binary = False)
        
        elif message['header'] == 'get genome':
            genome = self.clients[self.id]['db']['genomes'].find_one({'_id': message['genome_id']})
            answer = {'header': 'got genome'}
            answer['genome'] = genome
            self.write_message(answer, binary = False)
        
        elif message['header'] == 'get annotations per genome':
            annotations = []
            for annotation in self.clients[self.id]['db']['annotations'].find({'genome': message['genome_id']+"@genomes"}):
                annotations.append(annotation) 
            answer = {'header': 'got annotations per genome'}
            answer['annotations'] = annotations
            answer['center'] = message['center']
            answer['annotation_id'] = message['annotation_id']
            self.write_message(answer, binary = False)
        
        elif message['header'] == 'genome browser dragged':
            current_ids = message['current_ids']
            genomic_range = message['genomic_range']
            for annotation in self.clients[self.id]['db']['annotations'].find({'genome': message['genome_id']+"@genomes"}):
                if not annotation['_id'] in current_ids and annotation['genomicPositions'][0] >= genomic_range[0] and annotation['genomicPositions'][0] <= genomic_range[1] or annotation['genomicPositions'][1] >= genomic_range[0] and annotation['genomicPositions'][1] <= genomic_range[1]:
                    answer = {'header': 'got new annotation to display'}
                    answer['annotation'] = annotation
                    self.write_message(answer, binary = False)
        
    def on_close(self):
        self.clients.pop(self.id, None)
        print "Client disconnected"

class Application(tornado.web.Application):
    def __init__(self):

        handlers = [
            (r'/', IndexHandler),
            (r'/websocket', WebSocketHandler)
        ]

        settings = {
            'template_path': 'templates',
            'static_path': 'static'
        }

        tornado.web.Application.__init__(self, handlers, **settings)

if __name__ == "__main__":
    mongodb_host = "localhost"
    mongodb_port = 27017

    if "-mh" in sys.argv:
        mongodb_host = sys.argv[sys.argv.index("-mh")+1]
    if "-mp" in sys.argv:
        mongodb_port = int(sys.argv[sys.argv.index("-mp")+1])

    try :
        mongo_client = MongoClient(mongodb_host, mongodb_port)
    except Exception, e:
        print 'Cannot connect any Mongodb instance hosted at %s:%i'%(mongodb_host, mongodb_port)
        print 'Usage: ./server.py [-mh mongodb_host (default: localhost)] [-mp mongodb_port (default: 27017)]'
        sys.exit(-1)

    tornado.options.parse_command_line()
    app = Application()
    server = tornado.httpserver.HTTPServer(app)
    server.listen(8888)
    main_loop = tornado.ioloop.IOLoop.instance()
    main_loop.add_timeout(datetime.timedelta(seconds=5), webbrowser.open("http://localhost:8888"))
    main_loop.start()

