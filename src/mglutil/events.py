#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2011
#
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/events.py,v 1.4.2.1 2015/08/26 18:55:06 sanner Exp $
#
# $Id: events.py,v 1.4.2.1 2015/08/26 18:55:06 sanner Exp $
#
from time import time
import warnings

class Event:
    """
    Base class for events.
    """

    def __init__(self, *args, **kw):
        """  """
        self.timestamp = time()
        self.args = args
        self.kw = kw


    def __getattr__(self, name):
        return self.kw[name]

    
## helper functions to register and unregister listeners and dispatch events
def registerListener(eventClass, function, eventManager=None):
    assert issubclass(eventClass, Event)
    if eventManager:
        assert isinstance(eventManager, EventHandler)
    else:
        eventManager = _DefaultEventsManager
    eventManager.registerListener(eventClass, function)
    

def unregisterListener(eventClass, function, eventManager=None):
    assert issubclass(eventClass, Event)
    if eventManager:
        assert isinstance(eventManager, EventHandler)
    else:
        eventManager = _DefaultEventsManager
    eventManager.unregisterListener(eventClass, function)

    
def dispatchEvent(event, eventManager=None):
    assert isinstance(event, Event)
    if eventManager:
        assert isinstance(eventManager, EventHandler)
    else:
        eventManager = _DefaultEventsManager
    eventManager.dispatchEvent(event)
    

class EventHandler:
    """This class adds methods for registening functions called
listeners to be called upon a particular Event.
"""

    def __init__(self):
        self.eventListeners = {}


    def registerListener(self, event, function):
        """
        registers a function to be called for a given event.
        event has to be a class subclassing Event

        None <- registerListener(event, function)

        arguments:
            event: event class
            function: callable object that will be called with the
                      event instance as an argument.
        """
        assert issubclass(event, Event)
        assert callable(function)

        if event not in self.eventListeners:
            self.eventListeners[event] = [function]
        else:
            if function in self.eventListeners[event]:
                warnings.warn('function %s already registered for event %s'%(
                    function,event))
            else:
                self.eventListeners[event].append(function)


    def unregisterListener(self, event, function):
        """
        unregisters a function to be called for a given event.
        event has to be a class subclassing Event

        None <- unregisterListener(event, function)

        arguments:
            event: event class
            function: callable object that will be called with the
                      event instance as an argument.
        """
        if event in self.eventListeners:
            if function in self.eventListeners[event]:
                self.eventListeners[event].remove(function)

                
    def dispatchEvent(self, event):
        """call all registered listeners for this event type.
arguments:
    event: instance of an event
"""
        assert isinstance(event, Event)
        if event.__class__ in self.eventListeners:
            for func in self.eventListeners[event.__class__]:
                func(event)

try:
    _DefaultEventsManager
except NameError:
    _DefaultEventsManager = EventHandler()
    
## TODO add unregisterListener method
                
