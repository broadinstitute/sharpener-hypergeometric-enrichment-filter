import connexion
import six

from swagger_server.models.gene_info import GeneInfo  # noqa: E501
from swagger_server.models.transformer_info import TransformerInfo  # noqa: E501
from swagger_server.models.transformer_query import TransformerQuery  # noqa: E501
from swagger_server import util
from swagger_server.models.parameter import Parameter
from swagger_server.models.gene_info import GeneInfoIdentifiers
from swagger_server.models.attribute import Attribute

import sys
import json
import requests

class Transformer:

    def __init__(self, variables):
        with open("transformer_info.json",'r') as f:
            self.info = TransformerInfo.from_dict(json.loads(f.read()))
            self.variables = variables
            self.parameters = dict(zip(variables, self.info.parameters))


    def transform(self, query):
        query_controls = {control.name: control.value for control in query.controls}
        controls = {}
        for variable, parameter in self.parameters.items():
            if parameter.name in query_controls:
                controls[variable] = Transformer.get_control(query_controls[parameter.name], parameter)
            else:
                msg = "required parameter '{}' not specified".format(parameter.name)
                return ({ "status": 400, "title": "Bad Request", "detail": msg, "type": "about:blank" }, 400 )

        if self.info.function == 'producer':
            return self.produce(controls)
        if self.info.function == 'expander':
            return self.expand(query.genes, controls)
        if self.info.function == 'filter':
            return self.filter(query.genes, controls)

        return ({ "status": 500, "title": "Internal Server Error", "detail": self.info.function+" not implemented", "type": "about:blank" }, 500 )

    def produce(self, controls):
        return ({ "status": 500, "title": "Internal Server Error", "detail": "Producer not implemented", "type": "about:blank" }, 500 )


    def expand(self, query_genes, controls):
        return ({ "status": 500, "title": "Internal Server Error", "detail": "Expander not implemented", "type": "about:blank" }, 500 )


    def filter(self, query_genes, controls):
        return ({ "status": 500, "title": "Internal Server Error", "detail": "Filter not implemented", "type": "about:blank" }, 500 )


    @staticmethod
    def get_control(value, parameter):
        if parameter.type == 'double':
            return float(value)
        elif parameter.type == 'Boolean':
            return bool(value)
        elif parameter.type == 'int':
            return int(value)
        else:
            return value


