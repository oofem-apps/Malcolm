class element:
    def __init__(self, globvar, nr, nodes, cs, tag = None):
        self.nr = nr
        self.nodes = nodes
        self.cs = cs
        self.tag = tag

        # check if the same cross-section
        for i_cs in globvar.cross_sections:
            if (i_cs.nr == cs):
                i_cs.elements.append(nr)

    
class truss_elem(element):
    def __init__(self, globvar, nr, nodes, cs, tag):
        super().__init__(globvar, nr, nodes, cs, tag)

        self.nr_nodes = 2
