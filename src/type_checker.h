#ifndef SIMIT_TYPE_CHECKER_H
#define SIMIT_TYPE_CHECKER_H

#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <exception>

#include "hir.h"
#include "hir_visitor.h"
#include "types.h"
#include "domain.h"
#include "program_context.h"
#include "error.h"
#include "ir.h"

namespace simit {
namespace hir {

class TypeChecker : public HIRVisitor {
public:
  TypeChecker(std::vector<ParseError> *errors) : 
    retField(ir::Field("", ir::Type())), 
    skipCheckDeclared(false),
    errors(errors) {}

  void typeCheck(Program::Ptr program) {
    program->accept(this);
  }

private:
  virtual void visit(RangeIndexSet::Ptr);
  virtual void visit(SetIndexSet::Ptr);
  virtual void visit(DynamicIndexSet::Ptr);
  virtual void visit(ElementType::Ptr);
  virtual void visit(Endpoint::Ptr);
  virtual void visit(SetType::Ptr);
  virtual void visit(TupleType::Ptr);
  virtual void visit(ScalarTensorType::Ptr);
  virtual void visit(NonScalarTensorType::Ptr);
  virtual void visit(Field::Ptr);
  virtual void visit(ElementTypeDecl::Ptr);
  virtual void visit(IdentDecl::Ptr);
  virtual void visit(ExternDecl::Ptr);
  virtual void visit(FuncDecl::Ptr);
  virtual void visit(VarDecl::Ptr);
  virtual void visit(ConstDecl::Ptr);
  virtual void visit(WhileStmt::Ptr);
  virtual void visit(IfStmt::Ptr);
  virtual void visit(IndexSetDomain::Ptr);
  virtual void visit(RangeDomain::Ptr);
  virtual void visit(ForStmt::Ptr);
  virtual void visit(PrintStmt::Ptr);
  virtual void visit(AssignStmt::Ptr);
  virtual void visit(ExprParam::Ptr);
  virtual void visit(MapExpr::Ptr);
  virtual void visit(OrExpr::Ptr);
  virtual void visit(AndExpr::Ptr);
  virtual void visit(XorExpr::Ptr);
  virtual void visit(EqExpr::Ptr);
  virtual void visit(NotExpr::Ptr);
  virtual void visit(AddExpr::Ptr);
  virtual void visit(SubExpr::Ptr);
  virtual void visit(MulExpr::Ptr);
  virtual void visit(DivExpr::Ptr);
  virtual void visit(ElwiseMulExpr::Ptr);
  virtual void visit(ElwiseDivExpr::Ptr);
  virtual void visit(NegExpr::Ptr);
  virtual void visit(ExpExpr::Ptr);
  virtual void visit(TransposeExpr::Ptr);
  virtual void visit(CallExpr::Ptr);
  virtual void visit(TensorReadExpr::Ptr);
  virtual void visit(FieldReadExpr::Ptr);
  virtual void visit(VarExpr::Ptr);
  virtual void visit(IntLiteral::Ptr);
  virtual void visit(FloatLiteral::Ptr);
  virtual void visit(BoolLiteral::Ptr);
  virtual void visit(DenseIntVector::Ptr);
  virtual void visit(DenseFloatVector::Ptr);
  virtual void visit(DenseNDTensor::Ptr);
  virtual void visit(DenseTensorLiteral::Ptr);
  //virtual void visit(Test::Ptr);

  typedef std::vector<ir::Type> Type;
  typedef std::shared_ptr<Type> TypePtr;
  typedef std::shared_ptr<ir::IndexSet> IndexSetPtr;

  class DimError : public std::exception {
    const char *what() const noexcept {
      return "mismatched dimension sizes";
    }
  };
  class TypeError : public std::exception {
    const char *what() const noexcept {
      return "mismatched element types";
    }
  };

  struct TensorValues {
    enum class Type {UNKNOWN, INT, FLOAT};

    TensorValues() : dimSizes(1), type(Type::UNKNOWN) {};

    void addDimension() { dimSizes.push_back(1); }
    void addIntValues(const unsigned len) {
      switch (type) {
        case Type::FLOAT:
          throw TypeError();
          break;
        case Type::UNKNOWN:
          type = Type::INT;
          break;
        default:
          break;
      }
      dimSizes[dimSizes.size() - 1] += len;
    }
    void addFloatValues(const unsigned len) {
      switch (type) {
        case Type::INT:
          throw TypeError();
          break;
        case Type::UNKNOWN:
          type = Type::FLOAT;
          break;
        default:
          break;
      }
      dimSizes[dimSizes.size() - 1] += len;
    }
    void merge(const TensorValues &other) {
      if (type != other.type) {
        throw TypeError();
      }
      if (dimSizes.size() - 1 != other.dimSizes.size()) {
        throw DimError();
      }
      dimSizes[dimSizes.size() - 1]++;
    }

    std::vector<unsigned> dimSizes;
    Type                      type;
  };
 
  void typeCheckVarOrConstDecl(VarDecl::Ptr, const bool = false);
  void typeCheckBinaryElwise(BinaryExpr::Ptr);
  void typeCheckBinaryBoolean(BinaryExpr::Ptr);

  void markCheckWritable(HIRNode::Ptr node) {
    if (isa<VarExpr>(node)) {
      checkWritable = node;
    } else if (isa<TensorReadExpr>(node)) {
      markCheckWritable(to<TensorReadExpr>(node)->tensor);
    } else if (isa<FieldReadExpr>(node)) {
      markCheckWritable(to<FieldReadExpr>(node)->setOrElem);
    }
  }

  TypePtr inferType(Expr::Ptr ptr) {
    retType.reset();
    ptr->accept(this);
    const TypePtr ret = retType;
    retType.reset();
    return ret;
  }
  IndexSetPtr getIndexSet(IndexSet::Ptr ptr) {
    retIndexSet.reset();
    ptr->accept(this);
    const IndexSetPtr ret = retIndexSet;
    retIndexSet.reset();
    return ret;
  }
  ir::Expr getExpr(Endpoint::Ptr ptr) {
    retExpr = ir::Expr();
    ptr->accept(this);
    const ir::Expr ret = retExpr;
    retExpr = ir::Expr();
    return ret;
  }
  ir::Type getIRType(hir::Type::Ptr ptr) {
    retIRType = ir::Type();
    ptr->accept(this);
    const ir::Type ret = retIRType;
    retIRType = ir::Type();
    return ret;
  }
  ir::Field getField(Field::Ptr ptr) {
    retField = ir::Field("", ir::Type());
    ptr->accept(this);
    const ir::Field ret = retField;
    retField = ir::Field("", ir::Type());
    return ret;
  }
  ir::Var getVar(IdentDecl::Ptr ptr) {
    retVar = ir::Var();
    ptr->accept(this);
    const ir::Var ret = retVar;
    retVar = ir::Var();
    return ret;
  }
  TensorValues getTensorVals(DenseTensorElement::Ptr ptr) {
    retTensorVals = TensorValues();
    ptr->accept(this);
    const TensorValues ret = retTensorVals;
    retTensorVals = TensorValues();
    return ret;
  }

  bool compareTypes(const ir::Type &l, const ir::Type &r) {
    return (l.kind() == r.kind() && l == r);
  }
  std::string typeString(const ir::Type &type) {
    std::stringstream oss;
    oss << type;
    if (type.isTensor() && type.toTensor()->isColumnVector) {
      oss << "'";
    }
    return oss.str();
  }
  std::string typeString(const TypePtr &type) {
    if (type->size() == 0) {
      return "none";
    }

    std::stringstream oss;
    oss << "\'";
    if (type->size() > 1) {
      oss << "(";
    }
    
    bool printDelimiter = false;
    for (const auto compType : *type) {
      if (printDelimiter) {
        oss << ", ";
      }
      oss << typeString(compType);
      printDelimiter = true;
    }

    if (type->size() > 1) {
      oss << ")";
    }
    oss << "\'";
    return oss.str();
  }
  void reportError(const std::string msg, HIRNode::Ptr loc) {
    const auto err = ParseError(loc->getLineBegin(), loc->getColBegin(), 
                                loc->getLineEnd(), loc->getColEnd(), msg);
    errors->push_back(err);
  }
  void reportUndeclared(const std::string type, const std::string ident, 
                        HIRNode::Ptr loc) {
    std::stringstream errMsg;
    errMsg << "undeclared " << type << " \'" << ident << "\'";
    reportError(errMsg.str(), loc);
  }
  void reportMultipleDefs(const std::string type, const std::string ident, 
                          HIRNode::Ptr loc) {
    std::stringstream errMsg;
    errMsg << "multiple definitions of " << type << " \'" << ident << "\'";
    reportError(errMsg.str(), loc);
  }

  TypePtr retType;
  IndexSetPtr retIndexSet;
  ir::Expr retExpr;
  ir::Type retIRType;
  ir::Field retField;
  ir::Var retVar;
  TensorValues retTensorVals;
 
  bool skipCheckDeclared;
  HIRNode::Ptr checkWritable;
  internal::ProgramContext ctx;
  std::vector<ParseError> *errors;
};

}
}

#endif

