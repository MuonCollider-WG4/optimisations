import g4bl_interface

class CurrentBlockG4BL:
    def __init__(self, current_block):
        self.current_block = current_block
        self.element_list = []
    
    def make_g4bl_element(self):
        """
        self.zcentre = zcentre
        self.length = length
        self.rmin = rmin
        self.rmax = rmax
        self.period = period
        self.nrepeats = nrepeats
        self.nsheets = nsheets
        self.current_density = current_density
        self.coil_list = []
        """
        """
        self.z_position = 0.0
        self.inner_radius = 0.0
        self.outer_radius = 0.0
        self.length = 0.0
        self.current = 0.0
        self.name = "my_coil"
        """
        self.element_list = []
        coil_name = ""
        for i in [-self.nrepeats, self.nrepeats+1]:
            g4bl_solenoid = g4bl_interface.Solenoid()
            g4bl_solenoid.z_position = self.current_block.zcentre+i*self.current_block.period
            g4bl_solenoid.length = self.current_block.length
            g4bl_solenoid.rmin = self.current_block.inner_radius
            g4bl_solenoid.rmax = self.current_block.outer_radius
            g4bl_solenoid.nsheets = self.current_block.nsheets
            g4bl_solenoid.name = f"solenoid_{i}"
            g4bl_solenoid.current = self.current_block.current_density
            self.element_list.append(g4bl_solenoid)
        return self.element_list