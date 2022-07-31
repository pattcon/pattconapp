from django.urls import path
from . import views
from django.conf import settings
from django.conf.urls.static import static

urlpatterns = [path('', views.index), path('analise', views.analise),
path('motiflogo', views.motiflogo), path('fastaspecies', views.fastaspecies), path('fastagen', views.fastageneral),
 path('downfastsp', views.downfastsp), path('downfastgen', views.downfastgen), path('inicio', views.index), path('novapag', views.novapag),
path('report', views.report), path('motbysp', views.motbyspecies), path('tutorial', views.tutorial),
path('seqtest', views.downsequencetest),path('contact', views.contact), path('downfastgenTxt', views.downfastgenTxt)

               ]