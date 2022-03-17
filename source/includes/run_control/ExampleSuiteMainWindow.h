#ifndef CSMP_EXAMPLE_SUITE_MAIN_WINDOW_H
#define CSMP_EXAMPLE_SUITE_MAIN_WINDOW_H

#include "ExampleSuiteMainWindow.h"
#include "ExampleSuite.h"

// TODO: these needs the QT library, perhaps combine with our GUI

namespace PL{

class ExampleSuiteMainWindow : public QMainWindow
{
  Q_OBJECT

 public:
  ExampleSuiteMainWindow( QWidget *parent = 0 );

  void SetSuite( csmp::ExampleSuite* );

 private slots:
  // auto connected UI signals
  void on_pushButtonRun_clicked();
  void on_actionCSMP_triggered();
  void on_actionExample_Suite_triggered();
  void on_tableWidgetCategories_itemClicked( QTableWidgetItem *item );
  void on_tableWidgetExamples_itemClicked( QTableWidgetItem *item );
  void RunLoadedExample();

 private:
  // load form
  void Initialize();
  void LoadCategoryMap();
  void LoadCategoryMapToForm();
  void LoadExampleMap( QTableWidgetItem* categoryItem );
  void LoadExampleMapToForm( QTableWidgetItem* categoryItem );
  void LoadExampleDetailsToForm( QTableWidgetItem* exampleItem );
  // example suite
  csmp::ExampleSuite* suite_;
  csmp::Example* example_;
  // categories & examples db
  std::map<QTableWidgetItem*,csmp::ExampleSuite::category_iterator> categoryMap_;
  std::map<QTableWidgetItem*,csmp::Example*> exampleMap_;

  // UI composite
  Ui::ExampleSuiteMainWindow ui_;

};

} // PL

#endif // CSMP_EXAMPLE_SUITE_MAIN_WINDOW_H
