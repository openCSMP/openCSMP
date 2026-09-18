// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ExampleSuiteMainWindow.h"

#ifdef QT_WIDGETS_LIB // this code needs the QT library

#include <QMessageBox>
#include <QString>

using namespace std;
using namespace csmp;

namespace PL {

ExampleSuiteMainWindow::ExampleSuiteMainWindow( QWidget *parent )
{
  // setting up UI from form
  ui_.setupUi( this );
}

void ExampleSuiteMainWindow::SetSuite( ExampleSuite* suite )
{
  suite_ = suite;
  Initialize();
}

void ExampleSuiteMainWindow::Initialize()
{
  LoadCategoryMapToForm();
  LoadExampleMapToForm( categoryMap_.begin()->first );
}

/// loads categories from suite into categoryMap_
void ExampleSuiteMainWindow::LoadCategoryMap()
{
  QString categoryName;
  categoryMap_.clear();
  for( ExampleSuite::category_iterator it = (*suite_).GetCategoriesBegin();
       it != (*suite_).GetCategoriesEnd(); ++it )
  {
    categoryName = (*it).first.c_str();
    QTableWidgetItem* newItem = new QTableWidgetItem( categoryName );
    newItem->setFlags( Qt::ItemIsEditable );
    categoryMap_.insert( make_pair( newItem, it ) );
  }
}

/// loads categoryMap_ to form table
void ExampleSuiteMainWindow::LoadCategoryMapToForm()
{
  LoadCategoryMap();
  ui_.tableWidgetCategories->clear();
  for( map<QTableWidgetItem*,ExampleSuite::category_iterator>::iterator it = categoryMap_.begin();
       it != categoryMap_.end(); ++it )
  {
    ui_.tableWidgetCategories->insertRow( 0 );
    ui_.tableWidgetCategories->setItem( 0,0,(*it).first );
  }
}

/// loads examples of category(by table item) into exampleMap_
void ExampleSuiteMainWindow::LoadExampleMap( QTableWidgetItem* categoryItem )
{
  QString exampleName;
  exampleMap_.clear();
  vector<Example*>* exampleVector = categoryMap_.find( categoryItem )->second->second;
  for( vector<Example*>::iterator it = exampleVector->begin();
       it != exampleVector->end(); ++it )
  {
    exampleName = (*suite_).StripExampleTitle( (*it)->GetTitle() ).c_str();
    QTableWidgetItem* newItem = new QTableWidgetItem( exampleName );
    newItem->setFlags( Qt::ItemIsEditable );
    exampleMap_.insert( make_pair( newItem, (*it) ));
  }
}

/// loads exampleMap_ to form table
void ExampleSuiteMainWindow::LoadExampleMapToForm( QTableWidgetItem* categoryItem )
{
  LoadExampleMap( categoryItem );
  ui_.tableWidgetExamples->clear();
  ui_.tableWidgetExamples->setRowCount( 0 );
  for( map<QTableWidgetItem*,Example*>::iterator it = exampleMap_.begin();
       it != exampleMap_.end(); ++it )
  {
    ui_.tableWidgetExamples->insertRow( 0 );
    ui_.tableWidgetExamples->setItem( 0,0,(*it).first );
  } 
}

/// loads example details to form widget
void ExampleSuiteMainWindow::LoadExampleDetailsToForm( QTableWidgetItem* exampleItem )
{
  example_ = exampleMap_.find( exampleItem )->second;
  ui_.textEditExample->clear();
  QString cache;
  cache = (*suite_).StripExampleTitle( example_->GetTitle() ).c_str();
  ui_.textEditExample->append( cache );
  cache = "Authors: ";
  for( list<string>::const_iterator it = example_->GetAuthorsBegin(); it != example_->GetAuthorsEnd(); ++it )
  {
    cache += (*it).c_str();
    cache.append( " " );
  }
  ui_.textEditExample->append( cache );
  ui_.textEditExample->append( "Description: ");
  for( list<string>::const_iterator it = example_->GetDescriptionsBegin(); it != example_->GetDescriptionsEnd(); ++it )
  {
    cache = "";
    cache.append( "- " );
    cache += (*it).c_str();
    ui_.textEditExample->append( cache );
  }
  ui_.textEditExample->append( "Requirements: ");
  for( list<string>::const_iterator it = example_->GetRequirementsBegin(); it != example_->GetRequirementsEnd(); ++it )
  {
    cache = "";
    cache.append( "- " );
    cache += (*it).c_str();
    ui_.textEditExample->append( cache );
  }
}

/// runs loaded example
void ExampleSuiteMainWindow::RunLoadedExample()
{
  suite_->OstreamDoubleUnderlined( suite_->StripExampleTitle( example_->GetTitle() ) + " started..." );
  example_->Run();
  suite_->OstreamDoubleUnderlined( suite_->StripExampleTitle( example_->GetTitle() ) + " done..." );
}

void ExampleSuiteMainWindow::on_pushButtonRun_clicked()
{
  RunLoadedExample();
}

void ExampleSuiteMainWindow::on_tableWidgetCategories_itemClicked( QTableWidgetItem *item )
{
  LoadExampleMapToForm( item );
}

void ExampleSuiteMainWindow::on_tableWidgetExamples_itemClicked( QTableWidgetItem *item )
{
  LoadExampleDetailsToForm( item );
}



void ExampleSuiteMainWindow::on_actionCSMP_triggered()
{
  QMessageBox::about(this, tr("About CSMP"),
                           tr("<p>CSMP++ 1.0</p>"
                              "<p>SKM</p>"));
}

void ExampleSuiteMainWindow::on_actionExample_Suite_triggered()
{
  QMessageBox::about(this, tr("About the Example Suite GUI"),
                           tr("<p>Philipp Lang, 2010</p>"));
}

} // PL

#endif // if the QT widgets library is defined
