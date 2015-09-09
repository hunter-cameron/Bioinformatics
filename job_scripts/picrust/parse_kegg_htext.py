
import argparse
import re

class KeggClass(object):

    def __init__(self, name, parent=None):
        self.name = name

        # class B will have parents, class A will not
        self.parent = parent

        self.children = []

    def __str__(self):
        return self.name

    def add_child(self, child):
        self.children.append(child)

    

class KeggPath(object):

    def __init__(self, path_id, path_name, parent):
        self.path_id = path_id
        self.path_name = path_name
        
        self.parent = parent

        self.nodes = []

    def __str__(self):
        return self.path_name + " " + self.path_id

    def add_node(self, node):
        self.nodes.append(node)


class KeggNode(object):

    def __init__(self, id, description, pathway):
        self.id = id
    
        self.description = description
        self.pathways = [pathway]

    def get_full_pathstr(self):

        path_arr = []
        for pathway in self.pathways:
            full_path = [pathway]

            while full_path[0].parent:
                full_path = [full_path[0].parent] + full_path

            path_arr.append(";".join([str(item) for item in full_path]))


        return "|".join(path_arr)

class HtextParser(object):
    
    def __init__(self, htext_file):
        self.htext_file = htext_file

        self.regexes = {
                # A is simple, they take the format 'A<b>name</b>'
                "A": re.compile("A\<b\>(?P<name>[^<]+)\<\/b\>"),

                # B is the same as A except with some spaces in front 'B  <b>name</b>'
                # some B are blank (used as empty lines)
                "B": re.compile("B *\<b\>(?P<name>[^<]+)\<\/b\>"),

                # C takes the form 'C    01200 Carbon metabolism [PATH:ko01200]'
                "C": re.compile("C[\s\d]*(?P<name>.+?) \[(?:PATH|BR):(?P<path_id>[^]]*)\]"),

                # D takes a bunch of different forms but that's ok, we will get it all
                "D": re.compile("D *(?P<ko_id>K\d+) *(?P<description>.+?)$")
                
                }

    def parse(self):
        """
        This is a really messy file so I need to do excellent exception handling

        The interesting lines all begin with A, B, C, or D so I'll look at that first.

        """

        # set the KO levels
        class_A = None
        class_B = None
        pathway = None

        nodes = {}
        with open(self.htext_file, 'r') as IN:
            for line in IN:
                # see if line is interesting by first character

                if line[0] == "A":
                    # this pattern is simple; just <b>name</b>

                    m = self.regexes["A"].match(line)
                    name = m.group("name")

                    class_A = KeggClass(name)

                elif line[0] == "B":
                    # skip blank lines -- 2 for B and \n character
                    if len(line) == 2:
                        continue
   
                    m = self.regexes["B"].match(line)
                    name = m.group("name")
                    
                    class_B = KeggClass(name, class_A)
                    class_A.add_child(class_B)

                elif line[0] == "C":
                   
                    m = self.regexes["C"].match(line)
                    path_name = m.group("name")
                    path_id = m.group("path_id")

                    pathway = KeggPath(path_id, path_name, parent=class_B)
                    class_B.add_child(pathway)

                elif line[0] == "D":
                    
                    m = self.regexes["D"].match(line)
                    ko_id = m.group("ko_id")
                    description = m.group("description")

                    try:
                        node = nodes[ko_id].pathways.append(pathway)
                    except KeyError:
                        nodes[ko_id] = KeggNode(ko_id, description, pathway=pathway)

        return nodes



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="makes a KO metadata file ascribing each KO to a pathway and description")
    parser.add_argument("-htext", help="the 'htext' file from their website")
    parser.add_argument("-out", help="output path", default="ko_metadata.tab")

    args = parser.parse_args()

    htp = HtextParser(args.htext)
    nodes = htp.parse()

    sorted_nodes = sorted(nodes.values(), key=lambda node: node.id)

    with open(args.out, 'w') as OUT:
        OUT.write("\t".join(["KO", "Description", "Pathways"]) + "\n")

        for node in sorted_nodes:
            OUT.write("\t".join([node.id, node.description, node.get_full_pathstr()]) + "\n")

