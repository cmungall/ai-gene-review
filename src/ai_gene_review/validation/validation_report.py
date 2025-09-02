"""Pydantic models for validation reports.

This module defines the data structures for validation results,
including support for different severity levels (hard/soft failures).
"""

from enum import Enum
from pathlib import Path
from typing import List, Optional, Dict, Any
from datetime import datetime
from pydantic import BaseModel, Field, ConfigDict


class ValidationSeverity(str, Enum):
    """Severity levels for validation issues."""

    ERROR = "error"  # Hard fail - must be fixed
    WARNING = "warning"  # Soft fail - should be fixed but not critical
    INFO = "info"  # Informational - suggestions for improvement


class ValidationIssue(BaseModel):
    """A single validation issue found during validation."""

    model_config = ConfigDict(use_enum_values=True)

    severity: ValidationSeverity = Field(description="Severity level of the issue")
    message: str = Field(description="Human-readable description of the issue")
    path: Optional[str] = Field(
        default=None, description="JSON path or field path where the issue occurred"
    )
    details: Optional[Dict[str, Any]] = Field(
        default=None, description="Additional details about the issue"
    )
    suggestion: Optional[str] = Field(
        default=None, description="Suggested fix for the issue"
    )

    def is_error(self) -> bool:
        """Check if this is a hard failure."""
        return self.severity == ValidationSeverity.ERROR

    def is_warning(self) -> bool:
        """Check if this is a soft failure."""
        return self.severity == ValidationSeverity.WARNING

    def __str__(self) -> str:
        """String representation of the issue."""
        parts = [f"[{self.severity.upper()}]"]
        if self.path:
            parts.append(f"at {self.path}")
        parts.append(f": {self.message}")
        if self.suggestion:
            parts.append(f" (Suggestion: {self.suggestion})")
        return " ".join(parts)


class ValidationReport(BaseModel):
    """Complete validation report for a file."""

    model_config = ConfigDict(use_enum_values=True)

    file_path: Optional[Path] = Field(
        default=None, description="Path to the validated file"
    )
    is_valid: bool = Field(
        description="Whether the file passes validation (no errors, warnings allowed)"
    )
    issues: List[ValidationIssue] = Field(
        default_factory=list, description="List of validation issues found"
    )
    timestamp: datetime = Field(
        default_factory=datetime.now, description="When the validation was performed"
    )
    schema_version: Optional[str] = Field(
        default=None, description="Version of the schema used for validation"
    )
    metadata: Dict[str, Any] = Field(
        default_factory=dict, description="Additional metadata about the validation"
    )

    @property
    def error_count(self) -> int:
        """Count of hard failures (errors)."""
        return sum(1 for issue in self.issues if issue.is_error())

    @property
    def warning_count(self) -> int:
        """Count of soft failures (warnings)."""
        return sum(1 for issue in self.issues if issue.is_warning())

    @property
    def info_count(self) -> int:
        """Count of informational messages."""
        return sum(
            1 for issue in self.issues if issue.severity == ValidationSeverity.INFO
        )

    @property
    def has_errors(self) -> bool:
        """Check if there are any hard failures."""
        return self.error_count > 0

    @property
    def has_warnings(self) -> bool:
        """Check if there are any soft failures."""
        return self.warning_count > 0

    def get_errors(self) -> List[ValidationIssue]:
        """Get all error-level issues."""
        return [issue for issue in self.issues if issue.is_error()]

    def get_warnings(self) -> List[ValidationIssue]:
        """Get all warning-level issues."""
        return [issue for issue in self.issues if issue.is_warning()]

    def get_info(self) -> List[ValidationIssue]:
        """Get all info-level issues."""
        return [
            issue for issue in self.issues if issue.severity == ValidationSeverity.INFO
        ]

    def add_issue(
        self,
        severity: ValidationSeverity,
        message: str,
        path: Optional[str] = None,
        details: Optional[Dict[str, Any]] = None,
        suggestion: Optional[str] = None,
    ) -> None:
        """Add a validation issue to the report."""
        issue = ValidationIssue(
            severity=severity,
            message=message,
            path=path,
            details=details,
            suggestion=suggestion,
        )
        self.issues.append(issue)
        # Update validity based on errors
        if severity == ValidationSeverity.ERROR:
            self.is_valid = False

    def summary(self) -> str:
        """Generate a human-readable summary of the validation report."""
        if self.file_path:
            header = f"Validation Report for {self.file_path}"
        else:
            header = "Validation Report"

        lines = [
            header,
            "=" * len(header),
            f"Status: {'✓ Valid' if self.is_valid else '✗ Invalid'}",
            f"Errors: {self.error_count}",
            f"Warnings: {self.warning_count}",
            f"Info: {self.info_count}",
        ]

        if self.schema_version:
            lines.append(f"Schema Version: {self.schema_version}")

        if self.issues:
            lines.append("\nIssues:")
            for issue in self.issues:
                lines.append(f"  {issue}")

        return "\n".join(lines)

    def to_dict(self) -> Dict[str, Any]:
        """Convert report to dictionary for serialization."""
        return {
            "file_path": str(self.file_path) if self.file_path else None,
            "is_valid": self.is_valid,
            "error_count": self.error_count,
            "warning_count": self.warning_count,
            "info_count": self.info_count,
            "issues": [
                {
                    "severity": issue.severity,
                    "message": issue.message,
                    "path": issue.path,
                    "details": issue.details,
                    "suggestion": issue.suggestion,
                }
                for issue in self.issues
            ],
            "timestamp": self.timestamp.isoformat(),
            "schema_version": self.schema_version,
            "metadata": self.metadata,
        }


class BatchValidationReport(BaseModel):
    """Report for validation of multiple files."""

    reports: List[ValidationReport] = Field(
        default_factory=list, description="Individual validation reports for each file"
    )
    timestamp: datetime = Field(
        default_factory=datetime.now,
        description="When the batch validation was performed",
    )

    @property
    def total_files(self) -> int:
        """Total number of files validated."""
        return len(self.reports)

    @property
    def valid_files(self) -> int:
        """Number of valid files (no errors)."""
        return sum(1 for report in self.reports if report.is_valid)

    @property
    def invalid_files(self) -> int:
        """Number of invalid files (has errors)."""
        return sum(1 for report in self.reports if not report.is_valid)

    @property
    def files_with_warnings(self) -> int:
        """Number of files with warnings."""
        return sum(1 for report in self.reports if report.has_warnings)

    @property
    def total_errors(self) -> int:
        """Total number of errors across all files."""
        return sum(report.error_count for report in self.reports)

    @property
    def total_warnings(self) -> int:
        """Total number of warnings across all files."""
        return sum(report.warning_count for report in self.reports)

    def get_invalid_reports(self) -> List[ValidationReport]:
        """Get reports for files that failed validation."""
        return [report for report in self.reports if not report.is_valid]

    def get_reports_with_warnings(self) -> List[ValidationReport]:
        """Get reports for files that have warnings."""
        return [report for report in self.reports if report.has_warnings]

    def summary(self) -> str:
        """Generate a summary of the batch validation."""
        lines = [
            "Batch Validation Summary",
            "========================",
            f"Total files: {self.total_files}",
            f"Valid: {self.valid_files}",
            f"Invalid: {self.invalid_files}",
            f"Files with warnings: {self.files_with_warnings}",
            f"Total errors: {self.total_errors}",
            f"Total warnings: {self.total_warnings}",
        ]

        if self.invalid_files > 0:
            lines.append("\nInvalid files:")
            for report in self.get_invalid_reports():
                if report.file_path:
                    lines.append(
                        f"  - {report.file_path} ({report.error_count} errors)"
                    )

        if self.files_with_warnings > 0:
            lines.append("\nFiles with warnings:")
            for report in self.get_reports_with_warnings():
                if report.file_path:
                    lines.append(
                        f"  - {report.file_path} ({report.warning_count} warnings)"
                    )

        return "\n".join(lines)
