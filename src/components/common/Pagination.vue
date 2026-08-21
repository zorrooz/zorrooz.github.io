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
  padding: var(--sp-12) 0 var(--sp-4);
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
  gap: var(--sp-2);
  justify-self: center;
}

.page-btn--nav {
  min-width: 34px;
  height: 34px;
  border: none;
  background: transparent;
  color: var(--fg-3);
  border-radius: 50%;
}

.page-btn--nav:hover {
  border-color: transparent;
  color: var(--primary);
}

.page-btn:disabled {
  opacity: 0.4;
  cursor: not-allowed;
}

.page-ellipsis {
  font-size: 14px;
  color: var(--fg-3);
  padding: 0 4px;
}
</style>
