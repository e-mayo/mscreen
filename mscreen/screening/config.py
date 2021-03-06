"""
Generate sample config file for each docking engine
Convert config file between docking engines
"""
import re

class ConfFile:
    def __ini__(self, file, backend):
        pass

class ConfReader:
    """
    Config file reader object. This is a factory-like class.
    It have 3 functions:
        register reader
        get_reader
        read_conf       
    You shoud create an ConfReader object then register each reader you want
    and then use read_conf to read the conf file. 
    """
    def __init__(self):
        self._readers = {}
    
    def register_reader(self, backend, reader):
        """
        This function add a reader to the _readers dictionary

        Parameters
        ----------
        backend : str
            Key of the reader.
        reader : Reader
            The conf file reader (Should inherit from Reader class).

        Returns
        -------
        None.

        """
        self._readers[backend] = reader

    def get_reader(self, backend):
        """
        Get reader object

        Parameters
        ----------
        backend : str
            Reader key.

        Raises
        ------
        ValueError
            If backend dont match any _readers keys.

        Returns
        -------
        reader : Reader
            Return the reader object.

        """
        reader = self._readers[backend]
        if not reader: raise ValueError(backend)
        return reader
    
    def read_conf(self,file_name,backend):
        """
        If you want to read a conf file use this function.
        It read the configuration file usinng the backend register 
        in _readers.
    

        Parameters
        ----------
        file_name : Path
            The conf file path.
        backend : str
            The backend to use. eg: 'vina'

        Returns
        -------
        Dict
            A dictionary containing the information of the conf file.
            
        Example:
        --------
        I: readeer.read_conf('conf.txt','vina')
        O: {'receptor': 'arg',
             'flex': 'arg',
             'ligand': 'arg',
             'center_x': 'arg',
             'center_y': 'arg'}

        """
        reader = self.get_reader(backend)
        conf = reader(file_name,backend)
        
        return conf.get_prop()
    
    
class ConfWriter:
    def __init__(self):
        self._writers = {}
    
    def register_writer(self, backend, writer):
        self._writers[backend] = writer

    def get_writer(self, backend):
        writer = self._writers[backend]
        if not writer: raise ValueError(backend)
        return writer
    
    def write_conf(self,keywords,backend,filename):
        writer = self.get_writer(backend)
        w = writer(keywords, backend, filename)
        w.write()
        return True
        
        



class Writer:
    def __init__(self,keywords, backend, filename):
        self.keywords = keywords        
        self.backend = backend
        self.filename = filename        
        self.lines = []
    
    def keyword_to_lines(self):
        pass 
    
    def write(self):    
        self.lines = self.keyword_to_lines()
        with open(self.filename,'w') as file:
            # for line in self.lines: 
                # file.write(line)
            file.writelines(self.lines)

class VinaWriter(Writer):
        
    def keyword_to_lines(self):        
        lines = []

        for key in self.keywords.keys():
            if self.keywords[key]:
                line = f'{key} = {self.keywords[key]}\n'
            else:
                line = f'{key}\n'
            lines.append(line)
        self.lines = lines
        return lines


class PlantsWriter(Writer):
    def keyword_to_lines(self):
        
        lines = []
        for key in self.keywords.keys():
            line = f'{key} {self.keywords[key]}\n'
            lines.append(line)
        self.lines = lines
        return lines


class LedockWriter(Writer):
    
    def keyword_to_lines(self):
        
        lines = []
        for key in self.keywords.keys():
            line = f'{key}\n'
            lines.append(line)
            if type(self.keywords[key]) == tuple or type(self.keywords[key]) == list:
                for value in self.keywords[key]:
                    line = f'{value}\n\n'
                    lines.append(line)
            elif type(self.keywords[key]) == str:
                line = f'{self.keywords[key]}\n\n'
                lines.append(line)
            else:
                print(type(self.keywords[key]),'#'*72)
                # print(self.keywords)
                raise ValueError(line)
        self.lines = lines
        return lines


class Reader:
    def __init__(self,file_name, backend):
        self.file_name = file_name        
        self.backend = backend
        self.text = self.read(file_name)
    
    def clean_text(self):
        pattern = r'#.*'
        p = re.compile(pattern,re.M)
        text = p.sub('\n',self.text)
        pattern = r'\n+'
        p = re.compile(pattern,re.M)
        return p.sub('\n',text)
        
    def read(self,file_name):
        with open(file_name) as file:
            text = file.read()
        return text   



class VinaReader(Reader):
    
    def get_prop(self):
        text = self.clean_text()
        self.ctext = text
        pattern = r'(\w*\s?)=(\s?-?\w*(?:\s|\n))'
        pattern = r'(\w*\s?)=\s?((?:\w+:\d+,?)+|(?:-?\w+.\w+)|\d)'
        pattern = r'^(\w+)\s*(?:=\s*([-\w._\\/-:,\s]+))?'
        pattern = r'^(\w+)\s*(?:=\s*([^\n]+))?'
        p = re.compile(pattern,re.M|re.I)
        props = p.findall(text)
        keywords = {}
        for prop in props:
            keywords[prop[0].strip()] = prop[1].strip()
        return keywords


  

class PlantsReader(Reader):
  
    def get_prop(self):
        text = self.clean_text()
        pattern = r'(^\w+\s+)((?:[\w\.]+\s?){,3})'
        pattern = r'^(\w+)\s*(?:\s*([\w._\\/-:,\s]+))?'
        pattern = r'^(\w+)\s*(?:\s*([^\n]+))?'
        p = re.compile(pattern,re.M)
        props = p.findall(text)
        keywords = {}
        for prop in props:
            keywords[prop[0].strip()] = prop[1].strip()
        return keywords



class LedockReader(Reader):    
        
        
    def get_prop(self):
        text = self.clean_text()
        pattern = r'receptor\n(.*)'
        p = re.compile(pattern,re.I|re.M)
        receptor = p.findall(text)[0].strip()
        
        pattern = r'Ligands list\n(.*)'
        p = re.compile(pattern,re.I|re.M)
        ligand_list = p.findall(text)[0].strip()
        
        pattern = r'Binding pocket\n(?P<x>.*)\n(?P<y>.*)\n(?P<z>.*)'
        p = re.compile(pattern,re.I|re.M)
        binding_pocket = p.findall(text)[0]
        
        pattern = r'RMSD\n(.*)'
        p = re.compile(pattern,re.I|re.M)
        rmsd = p.findall(text)[0].strip()
        
        pattern = r'Number of binding poses\n(.*)'
        p = re.compile(pattern,re.I|re.M)
        n_binding_poses = p.findall(text)[0].strip()
        
        keywords = {'Receptor': receptor,
                    'Ligands list': ligand_list,
                    'Binding pocket' :binding_pocket,
                    'RMSD': rmsd,
                    'Number of binding poses': n_binding_poses
                    }
        return keywords
    

    

reader = ConfReader()
reader.register_reader('vina', VinaReader)
reader.register_reader('ledock', LedockReader)
reader.register_reader('plants', PlantsReader)

writer = ConfWriter()
writer.register_writer('vina', VinaWriter)
writer.register_writer('ledock', LedockWriter)
writer.register_writer('plants', PlantsWriter)
# =============================================================================
# COVERT CONF FILE
# =============================================================================
class ConfConverter:
    
    @staticmethod
    def convert(keywords, from_backend,to_backend):
        
        if from_backend != 'vina' or to_backend != 'vina':
            converter = get_converter(from_backend,'vina')
            coords = converter(keywords)
            converter = get_converter('vina', to_backend)
            coords =  converter(keywords)
            return coords
        else:
            converter = get_converter(from_backend,'to_backend')
            coords = converter(keywords)
            return coords
    
    @staticmethod
    def get_converter(from_backend,to_backend):
        converter = {
            'vina':{'plants':coor_vina_to_plants,
                    'ledock':coor_vina_to_ledock
                    },
            'ledock':{'vina':coor_ledock_to_vina,
                      },
            'plants':{'vina':coor_plants_to_vina,
                      }
                    }
        return get_converter['from_backend']['to_backend']
        
            
            
    @staticmethod
    def coor_vina_to_plants(vina_keywords):
    
        x = float(vina_keywords['center_x'])
        y = float(vina_keywords['center_y'])
        z = float(vina_keywords['center_z'])
        
        dx = float(vina_keywords['size_x'])
        dy = float(vina_keywords['size_y'])
        dz = float(vina_keywords['size_z'])
        
        #r_2 = ((dz/2)**2 + ((dx/2)**2+(dy/2)**2))**0.5
        
        dxy = (dx**2 + dy**2)**0.5
        r = ((dz**2 + dxy**2)**0.5)/2
        
        plants_bindingsite = {}
        plants_bindingsite['bindingsite_center'] = f'{x} {y} {z}'
        plants_bindingsite['bindingsite_radius'] = r

        return plants_bindingsite
    
    @staticmethod
    def coor_vina_to_ledock(vina_keywords):
        
        x = float(vina_keywords['center_x'])
        y = float(vina_keywords['center_y'])
        z = float(vina_keywords['center_z'])
        
        dx = float(vina_keywords['size_x'])
        dy = float(vina_keywords['size_y'])
        dz = float(vina_keywords['size_z'])
        
        xmin = x - dx/2
        xmax = x + dx/2
    
        ymin = y - dy/2
        ymax = y + dy/2
        
        zmin = z - dz/2
        zmax = z + dz/2     
        
        ledock_bindingsite = {}   
    
        ledock_bindingsite['Binding pocket'] = f'{xmin} {xmax}',\
                                               f'{ymin} {ymax}',\
                                               f'{zmin} {zmax}'
       
        return ledock_bindingsite
    
    @staticmethod
    def coor_ledock_to_vina(ledock_keywords):
        
        x, y, z = [i.split(' ') for i in ledock_keywords['Binding pocket']]
        
        xmin, xmax = float(x[0]), float(x[1])
        ymin, ymax = float(y[0]), float(y[1])
        zmin, zmax = float(z[0]), float(z[1])
        
        dx = xmax - xmin
        dy = ymax - ymin
        dz = zmax - zmin
        
        vina_binding_site = {}
        
        vina_keywords['center_x'] = dx/2
        vina_keywords['center_y'] = dy/2
        vina_keywords['center_z'] = dz/2
        
        vina_keywords['size_x'] = dx
        vina_keywords['size_y'] = dy
        vina_keywords['size_z'] = dz
        
        return vina_binding_site
    
    @staticmethod   
    def coor_plants_to_vina(plants_keywords):
        
        x, y, z = [ float(i) for i in plants_keywords['bindingsite_center'].split(' ')]
        r = float(plants_keywords['bindingsite_radius'])
        
        vina_binding_site = {}
        
        vina_binding_site['center_x'] = x
        vina_binding_site['center_y'] = y
        vina_binding_site['center_z'] = z
        
        vina_binding_site['size_x'] = r*2
        vina_binding_site['size_y'] = r*2
        vina_binding_site['size_z'] = r*2
        
        return vina_binding_site

    @staticmethod
    def speed_coverter(keywords):
        exhaustiveness = ['speed1','speed2','speed4']
        search_speed = ['32','16','8']
        
        if 'exhaustiveness' in keywords:
            speed = keywords['exhaustiveness']
            if int(speed) < 16:
                speed = 'speed4'
            elif 16 <= int(speed) < 32:
                speed = 'speed2'
            elif 32 <= int(speed):
                speed = 'speed1'
            return {'search_speed':speed}
        elif 'search_speed' in keywords:
            if keywords['search_speed'] == 'speed1':
                speed = '32'
            elif keywords['search_speed'] == 'speed2':
                speed = '16'
            elif keywords['search_speed'] == 'speed4':
                speed = '8'
            return {'exhaustiveness':speed}
        else:
            return None
        
    #%%
       
    @staticmethod
    def get_keywords_converter(from_backend,to_backend):
        """uncompleted"""
        
        
        keywords_vina_ledock={
                            'receptor':'Receptor',
                            'num_modes':'Number of binding poses'                    
                            }
        keywords_vina_plants={
                            'ligand':'ligand_file',
                            'receptor':'protein_file',
                            'num_modes':'cluster_structures'                    
                            }
        keywords_plants_ledock={
                              'cluster_rmsd':'RMSD',
                              'Receptor':'protein_file',
                              'ligand_list':'Ligands list',
                              'cluster_structures':'Number of binding poses'
                              }
        
        converter = {
            'vina':{'plants':keywords_vina_plants,
                    'ledock':keywords_vina_ledock
                    },
            'ledock':{'vina':'ledock_to_vina',
                      'plants':'ledock_to_plants',
                      },
            'plants':{'vina':'plants_to_vina',
                      'ledock':keywords_plants_ledock
                      },
            'general':{'plants':'vina_to_plants',
                       'ledock':'general_to_ledock',
                       'vina':'general_to_vina'
                      }
                    }
    @staticmethod   
    def covert_keywords(keywords,from_backend,to_backend):
        """
        Convert keywords from one backend to other
        only suport the default keywords
        
        Returns
        -------
        dict. converted keywords
        """
        pass
        # converter = get_keywords_converter(from_backend,to_backend)
#%%


#%%
if __name__ =='__main__':
    # =============================================================================
    # Test LedockReader
    # =============================================================================
    print('#'*72)    
    print('#'*72)    
    s = LedockReader('config_ledock_sample.txt', 'ledock')
    # pattern = r'#.*\n'
    # p = re.compile(pattern,re.M)
    # text = p.sub('\n',s.text)
    # pattern = r'\n+'
    # p = re.compile(pattern,re.M)
    # text = p.sub('\n',text)
    # print(text)
    print(len(s.get_prop()))
    for  i in s.get_prop():
        print(i,'-',s.get_prop()[i])
    print('#'*72)    
    print('#'*72) 
    #%% =============================================================================
    # Test VinaReader
    # =============================================================================
    print('#'*72)    
    print('#'*72)    
    s = VinaReader('config_smina_sample.txt', 'plants')    
    # pattern = r'#.*\n'
    # p = re.compile(pattern,re.M)
    # text = p.sub('\n',s.text)
    # pattern = r'\n+'
    # p = re.compile(pattern,re.M)
    # text = p.sub('\n',text)
    # print(text)
    print(len(s.get_prop()))
    k = s.get_prop()
    for  i in k.keys():
        print(i,'-',k[i])
    k.keys()
    print(k.values())
    print('#'*72)    
    print('#'*72)      
    #%% =============================================================================
    # TEST PlantsReader
    # =============================================================================
    print('#'*72)    
    print('#'*72)    
    s = PlantsReader('config_plants_sample.txt', 'plants')
    # pattern = r'#.*\n'
    # p = re.compile(pattern,re.M)
    # text = p.sub('\n',s.text)
    # pattern = r'\n+'
    # p = re.compile(pattern,re.M)
    # text = p.sub('\n',text)
    # print(text)
    print(len(s.get_prop()))
    for  i in s.get_prop():
        print(i,'-',s.get_prop()[i])
    print('#'*72)    
    print('#'*72)     
    
    #%% =============================================================================
    # Test ConfReader
    # =============================================================================
    
    reader = ConfReader()
    reader.register_reader('vina', VinaReader)
    reader.register_reader('ledock', LedockReader)
    reader.register_reader('plants', PlantsReader)
    
 
    vina_keywords   = reader.read_conf('config_vina_sample.txt', 'vina')
    ledock_keywords = reader.read_conf('config_ledock_sample.txt', 'ledock')
    plants_keywords = reader.read_conf('config_plants_sample.txt', 'plants')
    
    keywords = plants_keywords
    print(len(keywords))
    for  i in keywords:
        print(i,'-',keywords[i])
    #%% =============================================================================
    # Test ConfWriter
    # =============================================================================
    print('# '+'='*70)
    print('#     Test Write')
    print('# '+'='*70)
    writer = ConfWriter()
    writer.register_writer('vina', VinaWriter)
    writer.register_writer('ledock', LedockWriter)
    writer.register_writer('plants', PlantsWriter)
    
    writer.write_conf(vina_keywords, 'vina', 'cvina_writer.txt')
    writer.write_conf(ledock_keywords, 'ledock', 'cledock_writer.txt')
    writer.write_conf(plants_keywords, 'plants', 'cplants_writer.txt')
    
    #%% =============================================================================
    # TEST VinaWriter    
    # =============================================================================
    w = VinaWriter(vina_keywords,'vina','conf_vina_writer.txt')
    w.write()
    #%% =============================================================================
    # TEST PlantsWriter
    # =============================================================================
    w = PlantsWriter(plants_keywords,'platns','conf_plants_writer.txt')
    w.write()
    #%% =============================================================================
    # TEST LedockWriter
    # =============================================================================
    w = LedockWriter(ledock_keywords,'ledock','conf_ledock_writer.txt')
    w.write()
    #%% =============================================================================
    # TEST speed_coverter    
    # =============================================================================
    ConfConverter.speed_coverter(vina_keywords)
    ConfConverter.speed_coverter(plants_keywords)
    #%%
    # =============================================================================
    # Test coordinates transform
    # =============================================================================
    vina_box = {'center_x':-0.2 ,
                'center_y':-0.4,
                'center_z':1.9 ,
                'size_x': 14.1,
                'size_y': 17.1,
                'size_z': 15.1}
    vina_box = {'center_x': 54.0,
                'center_y':-19.4,
                'center_z':38.4,
                'size_x':17.3,
                'size_y':13.0,
                'size_z':17.4}
    
    plants_box = ConfConverter.coor_vina_to_plants(vina_box)