from argparse import ArgumentParser
from math import radians
from os.path import abspath, isfile

from networkx import Graph


class IEEECDFParser(object):
    def __init__(self, source_filename):
        # print abspath(cdf_input)
        if isfile(source_filename) is True:
            try:
                source_file = open(source_filename, 'r')
            except IOError:
                raise IOError('cannot open source IEEE common data format file')
        else:
            raise IOError('cannot open source IEEE common data format file: file does not exist or cannot accessed')

        self.source_file = source_file

        self.title_map = {
            "date": {"start": 1, "end": 9, "format_func": str},
            "originator_name": {"start": 10, "end": 30, "format_func": str},
            "mva_base": {"start": 31, "end": 37, "format_func": float},
            "year": {"start": 38, "end": 42, "format_func": int},
            "season": {"start": 43, "end": 43, "format_func": lambda x: "winter" if x == "W" else "summer"},
            "case_id": {"start": 45, "end": 73, "format_func": str}
        }

        self.bus_data_map = {
            "bus_num": {"start": 0, "end": 4, "format_func": int},
            # format outline says this should start at the 7th char, i.e., index 6
            "bus_name": {"start": 5, "end": 17, "format_func": str},
            "area_number": {"start": 18, "end": 20, "format_func": int},
            "loss_zone_number": {"start": 20, "end": 23, "format_func": int},
            "bus_type": {"start": 24, "end": 26, "format_func": int},
            "V": {"start": 27, "end": 33, "format_func": float},
            "theta": {"start": 33, "end": 40, "format_func": lambda x: radians(float(x))},
            "load_p": {"start": 40, "end": 49, "format_func": float},
            "load_q": {"start": 49, "end": 59, "format_func": float},
            "gen_p": {"start": 59, "end": 67, "format_func": float},
            "gen_q": {"start": 67, "end": 75, "format_func": float},
            "base_kv": {"start": 76, "end": 83, "format_func": float},
            "desired_volts": {"start": 84, "end": 90, "format_func": float},
            "max_limit": {"start": 90, "end": 98, "format_func": float},
            "min_limit": {"start": 98, "end": 106, "format_func": float},
            "shunt_g": {"start": 106, "end": 114, "format_func": float},
            "shunt_b": {"start": 114, "end": 121, "format_func": float},
            "remote_bus_num": {"start": 123, "end": 127, "format_func": int}
        }

        self.branch_data_map = {
            "tap_bus_num": {"start": 0, "end": 4, "format_func": int},
            "z_bus_num": {"start": 6, "end": 9, "format_func": int},
            "load_flow_area": {"start": 10, "end": 12, "format_func": int},
            # format outline says this should end at the 14th char, i.e., index 14
            "loss_zone": {"start": 12, "end": 15, "format_func": int},
            "circuit": {"start": 16, "end": 16, "format_func": int},
            "type": {"start": 18, "end": 18, "format_func": lambda x: "t_line" if int(x) == 0 else "fixed_tap" if int(
                x) == 1 else "unsupported"},
            "R": {"start": 19, "end": 29, "format_func": float},
            "X": {"start": 29, "end": 40, "format_func": float},
            "B": {"start": 40, "end": 50, "format_func": float},
            "line_mva_rating_1": {"start": 50, "end": 55, "format_func": int},
            "line_mva_rating_2": {"start": 56, "end": 61, "format_func": int},
            "line_mva_rating_3": {"start": 62, "end": 67, "format_func": int},
            "control_bus_num": {"start": 68, "end": 72, "format_func": int},
            "side": {"start": 73, "end": 73,
                     "format_func": lambda x: "term" if int(x) == 0 else "tap_bus" if int(x) == 1 else "z_bus"},
            "final_ratio": {"start": 76, "end": 82, "format_func": float},
            "final_angle": {"start": 83, "end": 90, "format_func": float},
            "min_tap": {"start": 90, "end": 97, "format_func": float},
            "max_tap": {"start": 97, "end": 104, "format_func": float},
            "step_size": {"start": 105, "end": 111, "format_func": float},
            "min_limit": {"start": 112, "end": 119, "format_func": float},
            "max_limit": {"start": 119, "end": 126, "format_func": float}
        }

    def _parse_line(self, mapping, line):
        parsed_line = {}
        for data_name in mapping.keys():
            parse_info = mapping[data_name]
            if parse_info["start"] == parse_info["end"]:
                chunk = line[parse_info["start"]]
            else:
                chunk = line[parse_info["start"]:parse_info["end"]]
            try:
                parsed_line[data_name] = parse_info["format_func"](chunk)
            except:
                raise Exception('could not parse line')

        return parsed_line

    def _parse_title_line(self, line):
        return self._parse_line(self.title_map, line)

    def _parse_bus_data_line(self, line):
        return self._parse_line(self.bus_data_map, line)

    def _parse_branch_data_line(self, line):
        return self._parse_line(self.branch_data_map, line)

    def generate_networkx_graph(self):
        # variables to store entire sections
        title_data = None
        bus_data_start_line = None
        bus_data_end_line = None
        branch_data_start_line = None
        branch_data_end_line = None
        Gn = Graph()
        for line_number, line in enumerate(self.source_file):
            line = line.rstrip()
            if line == "":
                continue
            if title_data is None:
                title_data = self._parse_title_line(line)
                for prop in title_data.keys():
                    Gn.graph[prop] = title_data[prop]
            elif bus_data_start_line is None:
                if line[0:16] == "BUS DATA FOLLOWS":
                    bus_data_start_line = line_number
            elif bus_data_start_line is not None and bus_data_end_line is None:
                if line[0:4] == "-999":
                    bus_data_end_line = line_number
                else:
                    bus_data = self._parse_bus_data_line(line)
                    Gn.add_node(bus_data["bus_num"], attr_dict=bus_data)
            elif branch_data_start_line is None:
                # format outline says this should end at character 16, i.e., index 15
                if line[0:19] == "BRANCH DATA FOLLOWS":
                    branch_data_start_line = line_number
            elif branch_data_start_line is not None and branch_data_end_line is None:
                if line[0:4] == "-999":
                    branch_data_end_line = line_number
                else:
                    branch_data = self._parse_branch_data_line(line)
                    bus_a = branch_data["tap_bus_num"]
                    bus_b = branch_data["z_bus_num"]
                    Gn.add_edge(bus_a, bus_b, attr_dict=branch_data)
                    # assuming pi model
                    Gn.node[bus_a]["shunt_b"] += 0.5 * branch_data["B"]
                    Gn.node[bus_b]["shunt_b"] += 0.5 * branch_data["B"]

        self.title_data = title_data
        self.bus_data_start_line = bus_data_start_line
        self.bus_data_end_line = bus_data_end_line
        self.branch_data_start_line = branch_data_start_line
        self.branch_data_end_line = branch_data_end_line
        self.Gn = Gn
        return Gn

IEEE_CASE = IEEECDFParser('ieee118cdf.txt')
Gn = IEEE_CASE.generate_networkx_graph()
1