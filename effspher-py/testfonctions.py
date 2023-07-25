import unittest
from viriel import surd_viriel_finale


class MyTestCase(unittest.TestCase):

    def test_val_surdensite_viriel(self):
        self.assertEqual(int(surd_viriel_finale), 173)


if __name__ == '__main__':
    unittest.main()
