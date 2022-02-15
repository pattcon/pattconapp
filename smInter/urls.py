from django.urls import path
from . import views

urlpatterns = [path('', views.index), path('analise', views.analise),
path('motivores', views.motivores), path('fastaspecies', views.fastaspecies), path('fastagen', views.fastageneral),
 path('downfastsp', views.downfastsp), path('downfastgen', views.downfastgen), path('inicio', views.index), path('novapag', views.novapag),
path('report', views.report)

               ]