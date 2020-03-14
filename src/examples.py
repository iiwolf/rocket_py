from constants import *
from orbital_mechanics import *
import unittest

class TestExamples(unittest.TestCase):
    
    def test1(self):
        
        # givens
        m_initial = 1924654
        m_prop = 1660998
        m_payload = 29484

        # calcs
        # m_inert = m_initial - m_prop - m_payload
        # imf = intert_mass_fraction(m_inert,m_prop)
        # mr = m_initial / (m_initial - m_prop)

        f_inert,mr,m_inert = generic_mass_sizing(m_initial, m_prop, m_payload)

        # asserts
        self.assertAlmostEqual(m_inert, 234172, 1)
        self.assertAlmostEqual(f_inert, 0.124, 1)
        self.assertAlmostEqual(mr, 7.3,1)

        # for 50% payload mass, how do these change?
        # imf,mr,m_inert = generic_mass_sizing(m_initial, m_prop, 0.50 * m_payload)

        m_payload *= 0.5
        m_prop = mass_propellant(m_payload, f_inert, mr)
        m_inert = f_inert * m_prop / (1 - f_inert)
        m_initial = m_prop + m_payload + m_inert

        # asserts
        self.assertAlmostEqual(m_prop, 830499, 1)
        self.assertAlmostEqual(m_inert, 117086, 1)
        self.assertAlmostEqual(m_initial, 962327,1)

        

         


unittest.main()