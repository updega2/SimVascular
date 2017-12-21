#include "svSegSelectionWidget.h"
#include "ui_svSegSelectionWidget.h"
#include "ui_svLoftParamWidget.h"
#include "svLoftingUtils.h"

#include <QMessageBox>

svSegSelectionWidget::svSegSelectionWidget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::svSegSelectionWidget)
    , m_TableModel(NULL)
    , m_NumSampling(0)
    , m_ModelElement(NULL)
    , m_ModelType("")
    , m_LoftWidget(NULL)
{
    ui->setupUi(this);

    ui->tableView->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
//    ui->tableView->setSortingEnabled(true);

    m_NodeMenu = new QMenu(ui->tableView);

    QAction* useSelectedAction=m_NodeMenu->addAction("Use Selected");
    QAction* useAllAction=m_NodeMenu->addAction("Use All");
    QAction* notUseAction=m_NodeMenu->addAction("Not Use Selected");
    QAction* useNoneAction=m_NodeMenu->addAction("Use None");

    QObject::connect( useSelectedAction, SIGNAL( triggered(bool) ) , this, SLOT( UseSelected(bool) ) );
    QObject::connect( useAllAction, SIGNAL( triggered(bool) ) , this, SLOT( UseAll(bool) ) );
    QObject::connect( notUseAction, SIGNAL( triggered(bool) ) , this, SLOT( NotUseSelected(bool) ) );
    QObject::connect( useNoneAction, SIGNAL( triggered(bool) ) , this, SLOT( UseNone(bool) ) );

    QObject::connect( ui->tableView, SIGNAL(customContextMenuRequested(const QPoint&))
      , this, SLOT(TableViewContextMenuRequested(const QPoint&)) );

    connect(ui->buttonBox,SIGNAL(accepted()), this, SLOT(Confirm()));
    connect(ui->buttonBox,SIGNAL(rejected()), this, SLOT(Cancel()));

    m_LoftWidget=new svLoftParamWidget();
    m_LoftWidget->move(400,400);
    m_LoftWidget->hide();
    m_LoftWidget->setWindowFlags(Qt::WindowStaysOnTopHint);

    connect(m_LoftWidget->ui->btnOK, SIGNAL(clicked()), this, SLOT(OKLofting()) );
    connect(m_LoftWidget->ui->btnApply, SIGNAL(clicked()), this, SLOT(ApplyLofting()) );
    connect(m_LoftWidget->ui->btnClose, SIGNAL(clicked()), this, SLOT(HideLoftWidget()) );

    connect(ui->btnUniformParameters, SIGNAL(clicked()), this, SLOT(ShowLoftWidget()) );
    connect(ui->btnSelectAllSegs, SIGNAL(clicked()), this, SLOT( UseAll() ));
    connect(ui->btnDeselectAllSegs, SIGNAL(clicked()), this, SLOT( UseNone() ));

}

svSegSelectionWidget::~svSegSelectionWidget()
{
    delete ui;

    if(m_LoftWidget)
        delete m_LoftWidget;
}

void svSegSelectionWidget::SetTableView(std::vector<mitk::DataNode::Pointer> segNodes, svModelElement* modelElement, std::string type)
{
    m_ModelElement=modelElement;
    m_ModelType=type;

    int segNum=segNodes.size();

    m_TableModel = new QStandardItemModel(segNum,2,this);

    for(int row = 0; row < segNum; row++)
    {
        for(int col = 0; col < 2; col++)
        {
            if(col==0)
            {
                QStandardItem* item= new QStandardItem(QString::fromStdString(segNodes[row]->GetName()));
                item->setEditable(false);
                m_TableModel->setItem(row,col,item);
            }
            else if(col==1)
            {
                if(modelElement&&modelElement->HasSeg(segNodes[row]->GetName()))
                {
                    QStandardItem* item= new QStandardItem(true);
                    item->setCheckable(true);
                    item->setCheckState(Qt::Checked);
                    m_TableModel->setItem(row,col,item);
                }
                else
                {
                    QStandardItem* item= new QStandardItem(false);
                    item->setCheckable(true);
                    item->setCheckState(Qt::Unchecked);
                    m_TableModel->setItem(row,col,item);
                }

            }
        }
    }

    QStringList headers;
    headers << "Segmentation" << "Use";
    m_TableModel->setHorizontalHeaderLabels(headers);

    ui->tableView->setModel(m_TableModel);
    //ui->tableView->setColumnWidth(0,150);

    int numSampling=0;
    if(modelElement)
    {
        numSampling=modelElement->GetNumSampling();

        if(numSampling>0)
            ui->lineEditNumSampling->setText(QString::number(numSampling));
        else
            ui->lineEditNumSampling->setText("");
    }
    else
    {
        if(type=="PolyData")
            ui->lineEditNumSampling->setText("");
        else
            ui->lineEditNumSampling->setText("20");
    }

    if(modelElement)
       ui->checkboxUseUniform->setChecked(modelElement->IfUseUniform()!=0);
    else
       ui->checkboxUseUniform->setChecked(false);

    if(type=="Parasolid")
        ui->widget_4->hide();
    else
        ui->widget_4->show();
}

std::vector<std::string> svSegSelectionWidget::GetUsedSegNames()
{
    std::vector<std::string> segNames;
    if(m_TableModel==NULL)
        return segNames;

//    int rowCount=m_TableModel->rowCount(QModelIndex());
    int rowCount=m_TableModel->rowCount();
    for(int i=0;i<rowCount;i++)
    {
        QModelIndex index= m_TableModel->index(i,1, QModelIndex());
        if(index.data(Qt::CheckStateRole) == Qt::Checked){
            QModelIndex idx= m_TableModel->index(i,0, QModelIndex());
            segNames.push_back(idx.data().toString().toStdString());
        }
    }
    return segNames;
}

int svSegSelectionWidget::GetNumSampling()
{
    return m_NumSampling;
}

int svSegSelectionWidget::IfUseUniform()
{
    return ui->checkboxUseUniform->isChecked()?1:0;
}

svLoftingParam svSegSelectionWidget::GetLoftingParam()
{
    return m_Param;
}

void svSegSelectionWidget::Confirm()
{
    QString strNum=ui->lineEditNumSampling->text().trimmed();
//    m_UseUniform = ui->checkboxUseUniform->isChecked()?1:0;

    if(strNum=="")
    {
        if(m_ModelType!="PolyData")
        {
            QMessageBox::warning(this,"Value Mising","Pleases provide the number of sampling points.");
            return;
        }

        m_NumSampling=0;
    }
    else
    {
        bool ok;
        int num=strNum.toInt(&ok);
        if(ok)
        {
            if(num<1)
            {
                QMessageBox::warning(this,"Value Error","Pleases give a positive integer format if you want to provide the number of sampling points.");
                return;
            }
            m_NumSampling=num;
        }
        else
        {
            QMessageBox::warning(this,"Format Error","Pleases give a correct format if you want to provide the number of sampling points.");
            return;
        }
    }
    hide();
    emit accepted();
}

void svSegSelectionWidget::Cancel()
{
    hide();
}

void svSegSelectionWidget::TableViewContextMenuRequested( const QPoint & pos )
{
    m_NodeMenu->popup(QCursor::pos());
}

void svSegSelectionWidget::UseSelected(bool)
{
    if(m_TableModel==NULL)
        return;

    QModelIndexList indexesOfSelectedRows = ui->tableView->selectionModel()->selectedRows();
    if(indexesOfSelectedRows.size() < 1)
    {
      return;
    }

    for (QModelIndexList::iterator it = indexesOfSelectedRows.begin()
       ; it != indexesOfSelectedRows.end(); it++)
     {
       int row=(*it).row();

       QStandardItem* item= m_TableModel->item(row,1);
       item->setCheckState(Qt::Checked);
     }
}

void svSegSelectionWidget::UseAll(bool)
{
    if(m_TableModel==NULL)
        return;

    int rowCount=m_TableModel->rowCount();

    for (int i=0;i<rowCount;i++)
    {
        QStandardItem* item= m_TableModel->item(i,1);
        item->setCheckState(Qt::Checked);
    }
}


void svSegSelectionWidget::NotUseSelected(bool)
{
    if(m_TableModel==NULL)
        return;

    QModelIndexList indexesOfSelectedRows = ui->tableView->selectionModel()->selectedRows();
    if(indexesOfSelectedRows.size() < 1)
    {
      return;
    }

    for (QModelIndexList::iterator it = indexesOfSelectedRows.begin()
       ; it != indexesOfSelectedRows.end(); it++)
     {
       int row=(*it).row();

       QStandardItem* item= m_TableModel->item(row,1);
       item->setCheckState(Qt::Unchecked);
     }
}

void svSegSelectionWidget::UseNone(bool)
{
    if(m_TableModel==NULL)
        return;

    int rowCount=m_TableModel->rowCount();

    for (int i=0;i<rowCount;i++)
    {
        QStandardItem* item= m_TableModel->item(i,1);
        item->setCheckState(Qt::Unchecked);
    }
}

void svSegSelectionWidget::ShowLoftWidget()
{
    if(m_ModelElement && m_ModelElement->IfUseUniform() && m_ModelElement->GetLoftingParam())
        m_Param=*(m_ModelElement->GetLoftingParam());
    else
        svLoftingUtils::SetPreferencedValues(&m_Param);

    m_LoftWidget->UpdateGUI(&m_Param);

    m_LoftWidget->show();
}

void svSegSelectionWidget::OKLofting()
{
    m_LoftWidget->UpdateParam(&m_Param);
    m_LoftWidget->hide();
}

void svSegSelectionWidget::ApplyLofting()
{
    m_LoftWidget->UpdateParam(&m_Param);
}

void svSegSelectionWidget::HideLoftWidget()
{
    m_LoftWidget->hide();
}
