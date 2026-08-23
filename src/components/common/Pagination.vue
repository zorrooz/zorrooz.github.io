<template>
  <nav class="pagination" v-if="totalPages > 1" :aria-label="t('pagination')">
    <div class="pagination__side">
      <button
        v-if="currentPage > 1"
        class="page-btn page-btn--nav"
        @click="onPrev?.()"
        :aria-label="t('prevPage')"
      >
        <i class="fas fa-chevron-left"></i>
      </button>
    </div>

    <div class="pagination__pages">
      <button
        v-if="showFirstPage"
        class="page-btn"
        :class="{ 'page-btn--active': currentPage === 1 }"
        @click="onGoToPage?.(1)"
      >
        1
      </button>
      <span class="page-ellipsis" v-if="showFirstEllipsis">...</span>

      <button
        v-for="page in middlePages"
        :key="page"
        class="page-btn"
        :class="{ 'page-btn--active': currentPage === page }"
        @click="onGoToPage?.(page)"
      >
        {{ page }}
      </button>

      <span class="page-ellipsis" v-if="showLastEllipsis">...</span>
      <button
        v-if="showLastPage && totalPages > 1"
        class="page-btn"
        :class="{ 'page-btn--active': currentPage === totalPages }"
        @click="onGoToPage?.(totalPages)"
      >
        {{ totalPages }}
      </button>
    </div>

    <div class="pagination__side pagination__side--right">
      <button
        v-if="currentPage < totalPages"
        class="page-btn page-btn--nav"
        @click="onNext?.()"
        :aria-label="t('nextPage')"
      >
        <i class="fas fa-chevron-right"></i>
      </button>
    </div>
  </nav>
</template>

<script setup lang="ts">
import { useI18n } from 'vue-i18n'

defineOptions({ name: 'PostPagination' })

interface PaginationProps {
  currentPage: number
  totalPages: number
  middlePages: number[]
  showFirstPage: boolean
  showLastPage: boolean
  showFirstEllipsis: boolean
  showLastEllipsis: boolean
  // eslint-disable-next-line no-unused-vars -- 类型签名无函数体，参数名仅作文档
  onGoToPage?: (page: number) => void
  onPrev?: () => void
  onNext?: () => void
}

defineProps<PaginationProps>()

const { t } = useI18n()
</script>

<style scoped>
.pagination {
  display: grid;
  grid-template-columns: 1fr auto 1fr;
  align-items: center;
  gap: var(--sp-2);
  margin-top: var(--sp-12);
  padding: var(--sp-8) 0 var(--sp-4);
  border-top: 1px solid var(--line);
}

.pagination__side {
  justify-self: start;
}

.pagination__side--right {
  justify-self: end;
}

.pagination__pages {
  display: flex;
  align-items: center;
  gap: var(--sp-1);
  justify-self: center;
}

/* 页码：方形文字按钮（radius 8 按钮规范），默认纯文字，当前页 tint 底 + 主色 */
.page-btn {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  min-width: 34px;
  height: 34px;
  padding: 0 6px;
  border: none;
  background: transparent;
  border-radius: var(--radius-btn);
  font-size: 13px;
  font-weight: 500;
  font-variant-numeric: tabular-nums;
  color: var(--fg-2);
  cursor: pointer;
  transition:
    color var(--dur-fast) ease,
    background-color var(--dur-fast) ease,
    transform var(--dur-fast) ease;
}

.page-btn:hover {
  color: var(--primary);
}

.page-btn:active {
  transform: scale(0.94);
}

.page-btn:focus-visible {
  outline: none;
  box-shadow: 0 0 0 2px var(--ring);
}

/* 当前页：主色 + tint 底（全局「选中背景」语义）。当前位置不可点——
   hover / focus / active 全部钉死为同一组值，任何交互态零变化 */
.page-btn--active,
.page-btn--active:hover,
.page-btn--active:focus-visible,
.page-btn--active:active {
  color: var(--primary);
  background: var(--tint);
  font-weight: 600;
  cursor: default;
  transform: none;
  box-shadow: none;
}

/* 前后翻页：圆形图标按钮（图标操作按钮规范） */
.page-btn--nav {
  width: 34px;
  padding: 0;
  border-radius: 50%;
  color: var(--fg-3);
  font-size: 12px;
}

.page-btn--nav:hover {
  color: var(--primary);
}

.page-btn:disabled {
  opacity: 0.4;
  cursor: not-allowed;
}

.page-ellipsis {
  font-size: 13px;
  color: var(--fg-3);
  padding: 0 4px;
  font-variant-numeric: tabular-nums;
}
</style>
