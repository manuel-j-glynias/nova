from jinja2 import Environment, BaseLoader
from datetime import datetime
import os



def render_all_pages(patient,dir_path):
    pages = ['page_1.html','page_2.html','page_3.html','page_5.html','page_5.html','page_6.html','page_7.html']
    wd = os.getcwd()
    static = wd + '/static/'
    for page in pages:
        with open('templates/'+page, 'r') as myfile:
            data = myfile.read()
            template = Environment(loader=BaseLoader()).from_string(data)
            # template.globals['STATIC_PREFIX'] = '/Users/mglynias/Desktop/PycharmProjects/nova/static/'
            template.globals['STATIC_PREFIX'] = static
            # output.append(template.render(patient=patient))
            with open(dir_path + '/'+ page, "w") as file:
                file.write(template.render(patient=patient))




def render_html_report(patient):
    order_id = patient['fake_order_id']
    now = datetime.now()
    timeStr = now.strftime("%m_%d_%Y__%H_%M")
    dir_path = 'output/' + order_id + '_' + timeStr
    os.mkdir(dir_path)
    render_all_pages(patient,dir_path)
