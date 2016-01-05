class FrozenClass(object):
    """ A class which can be (un-)frozen. If the class is frozen, no attributes
        can be added to the class. """

    __isfrozen = True

    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("%r is a frozen class. %r is not a valid \
attribute." % (self, key))

        object.__setattr__(self, key, value)

    def _freeze(self):
        """ Disallows adding new attributes. """
        self.__isfrozen = True

    def _unfreeze(self):
        """ Allows adding new attributes. """
        self.__isfrozen = False

    def _convert_type(self, k):
        attr = getattr(self, k)

        if isinstance(attr, bool):
            return attr
        try:
            return float(attr)
        except:
            return str(attr)

    def __str__(self):
        attrs = dir(self)
        attrs_dict = {}

        for k in attrs:
            if k.startswith("_"):
                continue
            val = self._convert_type(k)
            attrs_dict[k] = val

        return yaml.dump(attrs_dict, default_flow_style=False)
