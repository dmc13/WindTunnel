class BoundaryConditionSet(list):

    def add_bc(self, function_name, 
               expression=None, facet_id=None, bctype="strong_dirichlet"):
        """ Valid choices for bctype: "weak_dirichlet", "strong_dirichlet",
            "flather", "free_slip"
        """
        if expression is None and bctype!="free_slip":
            raise TypeError("Boundary condition of type %s requires "
                            "expression argument." % bctype)
        if expression is not None and bctype=="free_slip":
            raise TypeError('Boundary condition of type "free_slip" '
                            'does not allow expression argument.')
        if facet_id is None:
            raise TypeError('facet_id argument to add_bc() method is '
                            'not optional')

	self.append((function_name, expression, facet_id, bctype))
