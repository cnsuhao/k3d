#python

import testing

print """<DartMeasurement name="Test_Float" type="numeric/float">0.12</DartMeasurement>"""
print """<DartMeasurement name="Test_Integer" type="numeric/integer">12</DartMeasurement>"""
print """<DartMeasurement name="Test_String" type="text/string">A simple message</DartMeasurement>"""
print """<DartMeasurement name="Test_Text" type="text/text">A longer message</DartMeasurement>"""
print """<DartMeasurement name="Test_HTML" type="text/html"><![CDATA[<table><tr><th>foo</th><th>bar</th></tr><tr><td>a</td><td>b</td></tr></table>]]></DartMeasurement>"""
print """<DartMeasurement name="Test_XML" type="text/xml"><![CDATA[<foo><bar a="b"><baz c="d"/></bar></foo>]]></DartMeasurement>"""
print """<DartMeasurementFile name="Test_JPEG" type="image/jpeg">""" + testing.source_path() + """/bitmaps/test_rgb_8.jpg</DartMeasurementFile>"""
print """<DartMeasurementFile name="Test_PNG" type="image/png">""" + testing.source_path() + """/bitmaps/test_rgb_8.png</DartMeasurementFile>"""
